
#include <PeleLM.H>

#ifdef SPRAY_PELE_LM
#include "SprayParticles.H"

using namespace amrex;

namespace {
bool virtual_particles_set = false;
//
// Containers for the real "active" Particles
//
SprayParticleContainer* SprayPC = nullptr;
//
// Container for temporary, virtual Particles
//
SprayParticleContainer* VirtPC = nullptr;
//
// Container for temporary, ghost Particles
//
SprayParticleContainer* GhostPC = nullptr;

SprayData sprayData;
// Indices for spray source MultiFab
int sprayMomSrcIndx = VELX;
int sprayRhoSrcIndx = DENSITY;
int spraySpecSrcIndx = FIRSTSPEC;
int sprayEngSrcIndx = FIRSTSPEC + SPRAY_FUEL_NUM;
SprayComps scomps;

void
RemoveParticlesOnExit()
{
  delete SprayPC;
  SprayPC = nullptr;
  delete GhostPC;
  GhostPC = nullptr;
  delete VirtPC;
  VirtPC = nullptr;
}
std::string init_file;
int init_function = 1;
int spray_verbose = 0;
Real spray_cfl = 0.5;
Real wall_temp = 300.;
int mass_trans = 1;
int mom_trans = 1;
} // namespace

int PeleLM::do_spray_particles = 1;
// momentum + density + fuel species + enthalpy
int PeleLM::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;

int PeleLM::write_spray_ascii_files = 0;
int PeleLM::plot_spray_src = 0;
//Vector<std::string> PeleLM::spray_fuel_names;
std::string PeleLM::spray_fuel_names[SPRAY_FUEL_NUM];

SprayParticleContainer*
PeleLM::theSprayPC()
{
  return SprayPC;
}

SprayParticleContainer*
PeleLM::theVirtPC()
{
  return VirtPC;
}

SprayParticleContainer*
PeleLM::theGhostPC()
{
  return GhostPC;
}

Real
PeleLM::estSprayDt()
{
  Real estdt = 1.0e200;
  if (do_spray_particles != 1 || theSprayPC() == nullptr) {
    return estdt;
  }
  BL_PROFILE("PeleLM::sprayEstTimeStep()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    Real estdt_lev = theSprayPC()->estTimestep(lev, spray_cfl);
    if (estdt_lev > 0) {
      estdt = amrex::min(estdt, estdt_lev);
    }
  }
  return estdt;
}

void
PeleLM::readSprayParameters()
{
  ParmParse pp("peleLM");

  pp.query("do_spray_particles", do_spray_particles);
  if (do_spray_particles != 1) {
    do_spray_particles = 0;
    return;
  }
  SprayParticleContainer::readSprayParams(
    spray_verbose, spray_cfl, wall_temp, mass_trans, mom_trans, write_spray_ascii_files,
    plot_spray_src, init_function, init_file, sprayData, spray_fuel_names);
}

void
PeleLM::sprayParticleSetup()
{
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    amrex::Abort("Cannot have more spray fuel species than fluid species");
  }
#if NUM_SPECIES > 1
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == spray_fuel_names[i]) {
        sprayData.indx[i] = ns;
      }
    }
    if (sprayData.indx[i] < 0) {
      amrex::Print() << "Fuel " << spray_fuel_names[i]
                     << " not found in species list" << std::endl;
      amrex::Abort();
    }
  }
#else
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    sprayData.indx[ns] = 0;
  }
#endif
  amrex::Vector<Real> fuelEnth(NUM_SPECIES);
  auto eos = pele::physics::PhysicsType::eos();
  eos.T2Hi(sprayData.ref_T, fuelEnth.data());
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    const int fspec = sprayData.indx[ns];
    sprayData.latent[ns] -= fuelEnth[fspec] * 1.E-4;
  }
  scomps.mass_tran = mass_trans;
  scomps.mom_tran = mom_trans;
  // Component indices for conservative variables
  scomps.rhoIndx = DENSITY;
  scomps.momIndx = VELX;
  scomps.engIndx = RHOH;
  scomps.utempIndx = TEMP;
  scomps.specIndx = FIRSTSPEC;
  // Component indices for spray source MultiFab
  scomps.rhoSrcIndx = sprayRhoSrcIndx;
  scomps.momSrcIndx = sprayMomSrcIndx;
  scomps.specSrcIndx = spraySpecSrcIndx;
  scomps.engSrcIndx = sprayEngSrcIndx;
}

void
PeleLM::setupVirtualParticles(const int level)
{
  BL_PROFILE("PeleLM::setupVirtualParticles()");
  if (theSprayPC() != nullptr && !virtual_particles_set) {
    if (level < this->finestLevel()) {
#ifdef USE_SPRAY_SOA
      SprayParticleContainer::ParticleTileType virts;
#else
      SprayParticleContainer::AoS virts;
#endif
      this->setupVirtualParticles(level+1);
      theVirtPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);

      theSprayPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleLM::removeVirtualParticles(const int level)
{
  BL_PROFILE("PeleLM::removeVirtualParticles()");
  if (VirtPC != nullptr) {
    VirtPC->RemoveParticlesAtLevel(level);
  }
  virtual_particles_set = false;
}

void
PeleLM::setupGhostParticles(const int ngrow, const int level)
{
  BL_PROFILE("PeleLM::setupGhostParticles()");
  AMREX_ASSERT(level < this->finestLevel());
  if (theSprayPC() != nullptr) {
#ifdef USE_SPRAY_SOA
    SprayParticleContainer::ParticleTileType ghosts;
#else
    SprayParticleContainer::AoS ghosts;
#endif
    theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
    theGhostPC()->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleLM::removeGhostParticles(const int level)
{
  BL_PROFILE("PeleLM::removeGhostParticles()");
  if (GhostPC != nullptr) {
    GhostPC->RemoveParticlesAtLevel(level);
  }
}

void
PeleLM::createSprayData()
{
  SprayPC = new SprayParticleContainer(
    this, &m_phys_bc, sprayData, scomps, wall_temp);
  theSprayPC()->SetVerbose(spray_verbose);
  VirtPC = new SprayParticleContainer(
    this, &m_phys_bc, sprayData, scomps, wall_temp);
  GhostPC = new SprayParticleContainer(
    this, &m_phys_bc, sprayData, scomps, wall_temp);
}

void
PeleLM::initSprays()
{
  BL_PROFILE("PeleLM::initSprays()");

  //
  // Make sure to call RemoveParticlesOnExit() on exit.
  //
  amrex::ExecOnFinalize(RemoveParticlesOnExit);

  if (do_spray_particles) {
    if (theSprayPC() == nullptr) {
      createSprayData();
    }

    if (!init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(init_file, NSR_SPR + NAR_SPR);
    } else if (init_function > 0) {
      ProbParm const* lprobparm = prob_parm;
      theSprayPC()->InitSprayParticles(*lprobparm);
    } else {
      Abort("Must initialize spray particles with particles.init_function or "
            "particles.init_file");
    }
    if (spray_verbose >= 1) {
      amrex::Print() << "Total number of initial particles "
                     << theSprayPC()->TotalNumberOfParticles(false, false)
                     << std::endl;
    }
  }
}

void
PeleLM::sprayRestart()
{
  if (do_spray_particles) {
    AMREX_ASSERT(SprayPC == nullptr);
    createSprayData();

    //
    // Make sure to call RemoveParticlesOnExit() on exit.
    //
    amrex::ExecOnFinalize(RemoveParticlesOnExit);
    {
      theSprayPC()->Restart(
        m_restart_chkfile, "particles", true);
      amrex::Gpu::Device::streamSynchronize();
    }
  }
}

void
PeleLM::sprayMKD(const Real time,
                 const Real dt)
{
  if (spray_verbose) {
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";
  }
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (spray_verbose > 1) {
      amrex::Print() << "sprayMKDLevel " << lev << std::endl;
    }
    int spray_n_grow = 1;
    int ghost_width = 0;
    if (lev < finest_level) {
      ghost_width = 1;
    }
    int tmp_src_width = 1;
    if (lev > 0) tmp_src_width = 2;
    MultiFab tmp_spray_source(grids[lev], dmap[lev], num_spray_src, tmp_src_width, MFInfo(), Factory(lev));
    tmp_spray_source.setVal(0.);
    sprayMKDLevel(lev, time, dt, ghost_width, spray_n_grow, tmp_src_width, tmp_spray_source);
    theSprayPC()->transferSource(tmp_src_width, lev, tmp_spray_source, *(m_extSource[lev]));
  }
}

void
PeleLM::sprayMKDLevel(
  const int level,
  const Real time,
  const Real dt,
  const int ghost_width,
  const int spray_n_grow,
  const int tmp_src_width,
  MultiFab& tmp_spray_source)
{
  //
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  //
  int finest_level = this->finestLevel();

  //
  // Setup the virtual particles that represent particles on finer levels
  //
  if (level < finest_level) {
    setupVirtualParticles(level);
  }

  //
  // Make a copy of the particles on this level into ghost particles
  // for the finer level
  //
  if (level < finest_level) {
    setupGhostParticles(ghost_width, level);
  }

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  auto const* ltransparm = PeleLM::trans_parms.device_trans_parm();
  auto ldata_p = getLevelDataPtr(level,AmrOldTime);
  amrex::MultiFab& state = ldata_p->state;
  // Do the valid particles themselves
  theSprayPC()->moveKickDrift(
    state, tmp_spray_source, level, dt, time,
    false, // not virtual particles
    false, // not ghost particles
    spray_n_grow, tmp_src_width,
    true, ltransparm); // Move the particles

  // Only need the coarsest virtual particles here.
  if (level < finest_level) {
    theVirtPC()->moveKickDrift(
      state, tmp_spray_source, level, dt, time, true, false, spray_n_grow,
      tmp_src_width, true, ltransparm);
  }

  // Miiiight need all Ghosts
  if (theGhostPC() != nullptr && level != 0) {
    theGhostPC()->moveKickDrift(
      state, tmp_spray_source, level, dt, time, false, true, spray_n_grow,
      tmp_src_width, true, ltransparm);
  }
}

void
PeleLM::sprayMK(const Real time,
                const Real dt)
{
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (spray_verbose > 1) {
      amrex::Print() << "sprayMKLevel " << lev << std::endl;
    }
    int spray_n_grow = 1;
    int tmp_src_width = 1;
    if (lev > 0) tmp_src_width = 2;
    MultiFab tmp_spray_source(grids[lev], dmap[lev], num_spray_src, tmp_src_width, MFInfo(), Factory(lev));
    tmp_spray_source.setVal(0.);
    sprayMKLevel(lev, time, dt, spray_n_grow, tmp_src_width, tmp_spray_source);
    theSprayPC()->transferSource(tmp_src_width, lev, tmp_spray_source, *(m_extSource[lev]));
  }
}

void
PeleLM::sprayMKLevel(
  const int level,
  const Real time,
  const Real dt,
  const int spray_n_grow,
  const int tmp_src_width,
  amrex::MultiFab& tmp_spray_source)
{
  auto ldata_p = getLevelDataPtr(level,AmrNewTime);
  amrex::MultiFab& state = ldata_p->state;
  auto const* ltransparm = PeleLM::trans_parms.device_trans_parm();
  theSprayPC()->moveKick(
    state, tmp_spray_source, level, dt, time, false, false, spray_n_grow,
    tmp_src_width, ltransparm);

  if (level < this->finestLevel() && theVirtPC() != nullptr) {
    theVirtPC()->moveKick(
      state, tmp_spray_source, level, dt, time, true, false, spray_n_grow,
      tmp_src_width, ltransparm);
  }
  // Ghost particles need to be kicked except during the final iteration.
  if (theGhostPC() != nullptr && level != 0) {
    theGhostPC()->moveKick(
      state, tmp_spray_source, level, dt, time, false, true, spray_n_grow,
      tmp_src_width, ltransparm);
  }
}

void
PeleLM::sprayRedistribute(int lbase)
{
  if (lbase > 0) return;
  BL_PROFILE("PeleLM::sprayRedistribute");
  if (theSprayPC()) {
     static Vector<BoxArray> ba;
     static Vector<DistributionMapping> dm;
     bool changed = false;
     if (ba.size() != finest_level + 1) {
        ba.resize(finest_level + 1);
        dm.resize(finest_level + 1);
        changed = true;
     } else {
        for (int lev = 0; lev <= finest_level && !changed; lev++) {
           if (ba[lev] != grids[lev]) {
              changed = true;
           }
           if (!changed && dm[lev] != dmap[lev]) {
              changed = true;
           }
        }
     }
     if (changed) {
       theSprayPC()->Redistribute();
     }
  }
}

void
PeleLM::sprayPostTimestep()
{
  if (do_spray_particles) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      //
      // Remove virtual particles at this level if we have any.
      //
      removeVirtualParticles(lev);

      //
      // Remove Ghost particles on the final iteration
      //
      removeGhostParticles(lev);
      int nstep = 0; // Unused
      BL_PROFILE_VAR("SprayParticles::injectParticles()", INJECT_SPRAY);
      ProbParm const* lprobparm = prob_parm;
      Real cur_time = m_t_new[lev];
      Real dt = m_dt;
      bool injectParts = theSprayPC()->
        injectParticles(cur_time, dt, nstep, lev, finest_level, *lprobparm);
      BL_PROFILE_VAR_STOP(INJECT_SPRAY);
    }
    //
    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level.
    //
    theSprayPC()->Redistribute();
  }
}
#endif
