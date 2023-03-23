
#include "SprayParticles.H"
#include <AMReX_Particles.H>
#include "pelelm_prob.H"

using namespace amrex;

int
interpolateInjectTime(
  const Real& time, const int nvals, const Real* inject_time)
{
  int i = 0;
  Real ctime = inject_time[i];
  while (ctime < time) {
    ctime = inject_time[++i];
  }
  return i;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int nstep,
                                        int lev,
                                        int finest_level,
                                        ProbParm const& prob_parm)
{
  if (lev != 0) {
      return false;
  }

  if ( amrex::ParallelDescriptor::MyProc() == 0) {

    Real Pi_six = M_PI / 6.;

    // Get fuel species physical data
    const int pstateVel = m_sprayIndx.pstateVel;
    const int pstateT = m_sprayIndx.pstateT;
    const int pstateDia = m_sprayIndx.pstateDia;
    const int pstateY = m_sprayIndx.pstateY;
    const SprayData* fdat = m_sprayData;
    const Real num_ppp = fdat->num_ppp;
    const Real rho_part = fdat->rho[0];

    // Check how much mass we need to inject
    // and act accordingly
    Real mass_flow_rate = prob_parm.spray_mass_flow_rate;
    Real injection_mass = prob_parm.floating_injection_mass + mass_flow_rate * dt;

    // The following ensure that we don't mess with the size distribution and flow
    // rate, mostly critical for small dts. If we weren't able to inject the last particle
    // the distribution size returned, we stored it. Let's check again if we have enough mass
    // to inject it, if not just add whatever mass we need to inject to the internal reservoir.

    // This can't currently work because prob_parm is const !
    //Real part_mass_min = num_ppp Pi_six * rho_part * std::pow(prob_parm.part_dia_min, 3);
    //if ( injection_mass < part_mass_min ) {
    //   prob_parm.floating_injection_mass += mass_flow_rate * dt;
    //   return false;
    //}
    //if ( prob_parm.part_dia_saved > 0.0 ) {
    //  Real part_mass = num_ppp * Pi_six * rho_part * std::pow(prob_parm.part_dia_saved[inj], 3);
    //  if ( injection_mass < part_mass ) { // We can't do it, just add to mass counter and get out
    //    prob_parm.floating_injection_mass[inj] += mass_flow_rate * dt;
    //    return false;
    //  }
    //}

    // Geometry data
    const Geometry& geom = this->m_gdb->Geom(lev);
    const auto plo = geom.ProbLoArray();
    const auto phi = geom.ProbHiArray();
    const auto dx = geom.CellSize();
    RealVect dom_len(AMREX_D_DECL(geom.ProbLength(0),
                                  geom.ProbLength(1),
                                  geom.ProbLength(2)));
    AMREX_D_TERM(const Real splitx = plo[0] + 0.5 * dom_len[0];,
                 const Real splity = plo[1] + 0.5 * dom_len[1];,
                 const Real splitz = plo[2] + 0.5 * dom_len[2]);

    // Create a distribution
    std::unique_ptr<DistBase> dropDist;
    std::string dist_type = "Uniform";
    dropDist = DistBase::create(dist_type);
    dropDist->init("dropDist");

    // Host container
    amrex::ParticleLocData pld;
    std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>> host_particles;

    // Inject mass until we have the desired amount
    amrex::Real total_mass = 0.;
    amrex::Real remaining_mass = injection_mass;
    while (total_mass < injection_mass) {

      // Same comment as above: this won't work as long as prob_parm is const.
      // Get a particle size
      //Real cur_dia = 0.0;
      //if ( prob_parm.part_dia_saved[inj] > 0.0 ) { // We have a stored particle, use it (already checked that we can use it)
      //   cur_dia = prob_parm.part_dia_saved[inj];
      //   prob_parm.part_dia_saved[inj] = -1.0;
      //} else {    // Draw a new particle size
      //   cur_dia = dropDist->get_dia();
      //   Real pmass = Pi_six * rho_part * std::pow(cur_dia, 3);
      //   if ( remaining_mass < pmass ) {  // Don't have enough mass left, store it for later
      //       prob_parm.part_dia_saved[inj] = cur_dia;
      //       prob_parm.floating_injection_mass[inj] += remaining_mass;
      //       break;
      //   }
      //}
      Real cur_dia = dropDist->get_dia();

      // Injecting on x-hi
      RealVect part_loc(AMREX_D_DECL(phi[0],
                                     splity + (2.0 * amrex::Random() - 1.0) * prob_parm.jetRadius,
                                     0.0));

      ParticleType p;
      p.id() = ParticleType::NextID();
      p.cpu() = ParallelDescriptor::MyProc();

      // Particle velocity with a bit of random
      // NOTE: -xvel since injection on x+ face
      AMREX_D_TERM(Real x_vel = prob_parm.spray_vel + (2.0 * amrex::Random() - 1.0) * 0.05 * prob_parm.spray_vel ;,
                   Real y_vel = (2.0 * amrex::Random() - 1.0) * 0.05 * prob_parm.spray_vel ;,
                   Real z_vel = (2.0 * amrex::Random() - 1.0) * 0.05 * prob_parm.spray_vel );

      RealVect part_vel(AMREX_D_DECL(-x_vel, y_vel, z_vel));

      AMREX_D_TERM(p.rdata(pstateVel) = part_vel[0];,
                   p.rdata(pstateVel + 1) = part_vel[1];,
                   p.rdata(pstateVel + 2) = part_vel[2];);

      // Add particles as if they have advanced some random portion of dt
      Real pmov = amrex::Random();
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        p.pos(dir) = part_loc[dir] + pmov * dt * part_vel[dir];
      }

      p.rdata(pstateT) = prob_parm.spray_temp;
      p.rdata(pstateDia) = cur_dia;
      for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
        p.rdata(pstateY + sp) = 1.0;
      }

      // Put particle in place
      bool where = Where(p, pld);
      if (!where) {
          amrex::Abort("Bad injection particle");
      }
      std::pair<int, int> ind(pld.m_grid, pld.m_tile);

      host_particles[ind].push_back(p);

      Real pmass = Pi_six * rho_part * std::pow(cur_dia, 3);
      total_mass += num_ppp * pmass;
      remaining_mass = injection_mass - total_mass;
    }

    // Move particles to level holder
    for (auto& kv : host_particles) {
      auto grid = kv.first.first;
      auto tile = kv.first.second;
      const auto& src_tile = kv.second;
      auto& dst_tile = GetParticles(lev)[std::make_pair(grid, tile)];
      auto old_size = dst_tile.GetArrayOfStructs().size();
      auto new_size = old_size + src_tile.size();
      dst_tile.resize(new_size);
      // Copy the AoS part of the host particles to the GPU
      amrex::Gpu::copy( amrex::Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
                        dst_tile.GetArrayOfStructs().begin() + old_size);
    }
  }

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(const bool init_parts,
                                           ProbParm const& prob_parm)
{
  // This ensures the initial time step size stays reasonable
  m_injectVel = prob_parm.spray_vel;
  // Start without any particles
  return;
}
