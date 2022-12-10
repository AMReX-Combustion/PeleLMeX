
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

Real
jetOverlapArea(
  const Real testdx,
  const RealVect xloJ,
  const RealVect xhiJ,
  const RealVect jet_cent,
  const Real jr2,
#if AMREX_SPACEDIM == 3
  Real& loy,
  Real& hiy,
#endif
  Real& lox,
  Real& hix)
{
  // This is how much we reduce x to test for the area of the jet
  // This helps if the jet size is smaller than cell size
  Real cur_jet_area = 0.;
  Real testdx2 = testdx * testdx;
  Real curx = xloJ[0];
  hix = xloJ[0];
  lox = xhiJ[0];
  // Loop over each cell and check how much overlap there is with the jet
#if AMREX_SPACEDIM == 3
  hiy = xloJ[1];
  loy = xhiJ[1];
  while (curx < xhiJ[0]) {
    Real cury = xloJ[1];
    while (cury < xhiJ[1]) {
      Real r2 = curx * curx + cury * cury;
      if (r2 <= jr2) {
        cur_jet_area += testdx2;
        lox = amrex::min(curx, lox);
        hix = amrex::max(curx, hix);
        loy = amrex::min(cury, loy);
        hiy = amrex::max(cury, hiy);
      }
      cury += testdx;
    }
    curx += testdx;
  }
#else
  while (curx < xhiJ[0]) {
    Real r2 = curx * curx;
    if (r2 <= jr2) {
      cur_jet_area += testdx;
      lox = amrex::min(curx, lox);
      hix = amrex::max(curx, hix);
    }
    curx += testdx;
  }
#endif
  return cur_jet_area;
}

bool
SprayParticleContainer::injectParticles(Real time,
                                        Real dt,
                                        int nstep,
                                        int lev,
                                        int finest_level,
                                        ProbParm const& prob_parm)
{
  if (lev != 0)
      return false;

  // Get fuel species physical data
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateY = m_sprayIndx.pstateY;
  const SprayData* fdat = m_sprayData;
  Real rho_part = fdat->rho[0];

  // Check how much mass we need to inject
  // and act accordingly
  Real mass_flow_rate = prob_parm.mass_flow_rate;
  Real injection_mass = prob_parm.floating_injection_mass + mass_flow_rate * dt;
  Real Pi_six = M_PI / 6.;
  Real part_mass_min = Pi_six * rho_part * std::pow(prob_parm.part_dia_min, 3);
  if ( injection_mass < part_mass_min ) {
     // Don't have enough mass to inject a single part -> add to available mass
     // TODO: MPI ?
//     prob_parm.floating_injection_mass += mass_flow_rate * dt;
     return false;
  }

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

  // Injection data
  Real jet_vel = prob_parm.jet_vel;
  // This absolutely must be included with any injection or insertion
  // function or significant issues will arise
  if (jet_vel * dt / dx[0] > 0.5) {
    Real max_vel = dx[0] * 0.5 / dt;
    if (ParallelDescriptor::IOProcessor()) {
      std::string warn_msg =
        "Injection velocity of " + std::to_string(jet_vel) +
        " is reduced to maximum " + std::to_string(max_vel);
      amrex::Warning(warn_msg);
    }
    m_injectVel = jet_vel;
    jet_vel = max_vel;
  }
  // Injection on the x+ domain face. Get Area/Length of injection region
  // TODO: 2D ! it's twice the jetRadius
  Real jetArea = 2.*prob_parm.jetRadius;

  // Particle data
  const Real num_ppp = fdat->num_ppp;
  Real part_temp = prob_parm.part_temp;
  Real part_dia = prob_parm.part_mean_dia;
  Real part_stdev = prob_parm.part_stdev_dia;
  Real stdsq = part_stdev * part_stdev;
  Real meansq = part_dia * part_dia;
  Real log_mean = 2. * std::log(part_dia) - 0.5 * std::log(stdsq + meansq);
  Real log_stdev = std::sqrt(amrex::max(-2. * std::log(part_dia) + std::log(stdsq + meansq), 0.));

  for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();

    // Physical coordinates of the current box
    const RealBox& Rbox = RealBox(bx, geom.CellSize(), geom.ProbLo());
    const Real* xloB = Rbox.lo();
    const Real* xhiB = Rbox.hi();

    // Only work if intersection for x+ domain face
    if (xhiB[0] == phi[0]) {
      Gpu::HostVector<ParticleType> host_particles;
#ifdef USE_SPRAY_SOA
      std::array<Gpu::HostVector<Real>, NAR_SPR> host_real_attribs;
#endif

      // Get intersection of the current box with injection area to get
      // how much fuel need to be added in this box
      RealVect cur_jet_cent {AMREX_D_DECL(phi[0],splity,splitz)};  // Center of x+ face
      RealVect xlo;                                                // low/high intersection box/jet
      RealVect xhi;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          xlo[dir] = amrex::max(xloB[dir], cur_jet_cent[dir] - prob_parm.jetRadius);
          xhi[dir] = amrex::min(xhiB[dir], cur_jet_cent[dir] + prob_parm.jetRadius);
      }
      Real curJetArea = xhi[1] - xlo[1];

      if (curJetArea > 0.) {
        Real jet_perc = curJetArea / jetArea;
        Real perc_mass = jet_perc * injection_mass;
//        prob_parm.floating_injection_mass = 0.0;
        Real total_mass = 0.;
        Print() << " In this box, " << jet_perc * 100 << " % of the total " << injection_mass << " injected mass this dt \n";
        while (total_mass < perc_mass) {
            RealVect part_loc(AMREX_D_DECL(phi[0],
                                           xlo[1] + amrex::Random() * curJetArea,
                                           0.0));

            ParticleType p;
            p.id() = ParticleType::NextID();
            p.cpu() = ParallelDescriptor::MyProc();

            // Particle velocity with a bit of random
            // NOTE: -xvel since injection on x+ face
            AMREX_D_TERM(Real x_vel = jet_vel + (2.0 * amrex::Random() - 1.0) * 0.05 * jet_vel ;,
                         Real y_vel = (2.0 * amrex::Random() - 1.0) * 0.05 * jet_vel ;,
                         Real z_vel = (2.0 * amrex::Random() - 1.0) * 0.05 * jet_vel );
            RealVect part_vel(AMREX_D_DECL(-x_vel, y_vel, z_vel));
#ifdef USE_SPRAY_SOA
            AMREX_D_TERM(host_real_attribs[pstateVel].push_back(part_vel[0]);,
                         host_real_attribs[pstateVel + 1].push_back(part_vel[1]);,
                         host_real_attribs[pstateVel + 2].push_back(part_vel[2]));
#else
            AMREX_D_TERM(p.rdata(pstateVel) = part_vel[0];,
                         p.rdata(pstateVel + 1) = part_vel[1];,
                         p.rdata(pstateVel + 2) = part_vel[2];);
#endif

            Real cur_dia = amrex::RandomNormal(log_mean, log_stdev);
            // Use a log normal distribution
            cur_dia = std::exp(cur_dia);

            // Add particles as if they have advanced some random portion of dt
            Real pmov = amrex::Random();
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              p.pos(dir) = part_loc[dir] + pmov * dt * part_vel[dir];
            }
#ifdef USE_SPRAY_SOA
            host_real_attribs[pstateT].push_back(part_temp);
            host_real_attribs[pstateDia].push_back(cur_dia);
            for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp)
              host_real_attribs[pstateY + sp].push_back(1.0);
#else
            p.rdata(pstateT) = part_temp;
            p.rdata(pstateDia) = cur_dia;
            for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp)
              p.rdata(pstateY + sp) = 1.0;
#endif
            host_particles.push_back(p);
            Real pmass = Pi_six * rho_part * std::pow(cur_dia, 3);
            total_mass += num_ppp * pmass;
        }
      }

      // Move particles to level holder
      if (host_particles.size() > 0) {
        auto& particle_tile =
          GetParticles(lev)[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
        auto old_size = particle_tile.GetArrayOfStructs().size();
        auto new_size = old_size + host_particles.size();
        particle_tile.resize(new_size);

        Gpu::copy(
          Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
          particle_tile.GetArrayOfStructs().begin() + old_size);
#ifdef USE_SPRAY_SOA
        for (int i = 0; i != NAR_SPR; ++i) {
          Gpu::copy(
            Gpu::hostToDevice, host_real_attribs[i].begin(),
            host_real_attribs[i].end(),
            particle_tile.GetStructOfArrays().GetRealData(i).begin() +
              old_size);
        }
#endif
      }
    }
  }

  // Redistribute is done outside of this function
  return true;
}

void
SprayParticleContainer::InitSprayParticles(ProbParm const& prob_parm)
{
  // This ensures the initial time step size stays reasonable
  m_injectVel = prob_parm.jet_vel;
  // Start without any particles
  return;
}
