
#include <PeleLM.H>

#ifdef SPRAY_PELE_LM
#include "SprayParticles.H"

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
Real particle_cfl = 0.5;
Real wall_temp = 300.;
int mass_trans = 1;
int mom_trans = 1;
} // namespace

int PeleLM::do_spray_particles = 1;
// momentum + density + fuel species + enthalpy
int PeleLM::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;

int PeleLM::write_spray_ascii_files = 0;
int PeleLM::plot_spray_src = 0;
Vector<std::string> PeleLM::spray_fuel_names;

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

#if 0
void
PeleLM::particleEstTimeStep(Real& est_dt)
{
  if (do_spray_particles != 1 || theSprayPC() == nullptr) {
    return;
  }
  BL_PROFILE("PeleLM::particleEstTimeStep()");
  Real est_dt_particle = theSprayPC()->estTimestep(level, particle_cfl);

  if (est_dt_particle > 0) {
    est_dt = amrex::min(est_dt, est_dt_particle);
  }

  if (verbose && ParallelDescriptor::IOProcessor()) {
    if (est_dt_particle > 0) {
      amrex::Print() << "...estdt from particles at level " << level << ": "
                     << est_dt_particle << '\n';
    } else {
      amrex::Print() << "...there are no particles at level " << level << '\n';
    }
  }
}
#endif
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
    spray_verbose, particle_cfl, wall_temp, mass_trans, mom_trans, write_spray_ascii_files,
    plot_spray_src, init_function, init_file, sprayData, spray_fuel_names);
}

// Define gas phase state MF for computing spray source terms
#if 0
void
PeleLM::defineSprayStateMF()
{
  int nGrowS = 4;
  Sborder.define(grids, dmap, NUM_STATE, nGrowS, amrex::MFInfo(), Factory());
}
#endif
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
#if 0
void
PeleLM::setupVirtualParticles(int level)
{
  BL_PROFILE("PeleLM::setupVirtualParticles()");
  if (theSprayPC() != nullptr && !virtual_particles_set) {
    if (level < this->finestLevel()) {
#ifdef USE_SPRAY_SOA
      SprayParticleContainer::ParticleTileType virts;
#else
      SprayParticleContainer::AoS virts;
#endif
      ((PeleLM*)&parent->getLevel(level + 1))->setupVirtualParticles();
      theVirtPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);

      theSprayPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleLM::removeVirtualParticles()
{
  BL_PROFILE("PeleLM::removeVirtualParticles()");
  if (VirtPC != nullptr) {
    VirtPC->RemoveParticlesAtLevel(level);
  }
  virtual_particles_set = false;
}

void
PeleLM::setupGhostParticles(int ngrow, const int level)
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
#endif
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
    if (spray_verbose > 1) {
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
#if 0
void
PeleLM::particleMKD(
  const int level,
  const Real time,
  const Real dt,
  const int ghost_width,
  const int spray_n_grow,
  const int tmp_src_width,
  const int where_width,
  amrex::MultiFab& tmp_spray_source)
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
    setupVirtualParticles();
  }

  //
  // Make a copy of the particles on this level into ghost particles
  // for the finer level
  //
  if (level < finest_level) {
    setupGhostParticles(ghost_width);
  }

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  if (spray_verbose) {
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";
  }

  // Do the valid particles themselves
  theSprayPC()->moveKickDrift(
    Sborder, tmp_spray_source, level, dt, time,
    false, // not virtual particles
    false, // not ghost particles
    spray_n_grow, tmp_src_width,
    true, // Move the particles
    where_width);
  // Only need the coarsest virtual particles here.
  if (level < finest_level) {
    theVirtPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, true, false, spray_n_grow,
      tmp_src_width, true, where_width);
  }

  // Miiiight need all Ghosts
  if (theGhostPC() != nullptr && level != 0) {
    theGhostPC()->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, false, true, spray_n_grow,
      tmp_src_width, true, where_width);
  }
  MultiFab& spraydot = get_new_data(spraydot_Type);
  spraydot.setVal(0.);
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, spraydot);
}

void
PeleLM::particleMK(
  const Real time,
  const Real dt,
  const int spray_n_grow,
  const int tmp_src_width,
  amrex::MultiFab& tmp_spray_source)
{
  theSprayPC()->moveKick(
    Sborder, tmp_spray_source, level, dt, time, false, false, spray_n_grow,
    tmp_src_width);

  if (level < this->finestLevel() && theVirtPC() != nullptr) {
    theVirtPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, true, false, spray_n_grow,
      tmp_src_width);
  }
  // Ghost particles need to be kicked except during the final iteration.
  if (theGhostPC() != nullptr && level != 0) {
    theGhostPC()->moveKick(
      Sborder, tmp_spray_source, level, dt, time, false, true, spray_n_grow,
      tmp_src_width);
  }
  MultiFab& spraydot = get_new_data(spraydot_Type);
  spraydot.setVal(0.);
  theSprayPC()->transferSource(
    tmp_src_width, level, tmp_spray_source, spraydot);
}
#endif

void
PeleLM::sprayRedistribute(int lbase, bool init_part)
{
  if (do_spray_particles != 1) {
    return;
  }
  BL_PROFILE("PeleLM::sprayRedistribute()");
  int flev = this->finestLevel();
  if (theSprayPC()) {
    //
    // If we are calling with init_part = true, then we want to force the
    // redistribute without checking whether the grids have changed.
    //
    if (init_part) {
      theSprayPC()->Redistribute(lbase);
      return;
    }

    //
    // These are usually the BoxArray and DMap from the last regridding.
    //
    static Vector<BoxArray> ba;
    static Vector<DistributionMapping> dm;

    bool changed = false;

    if (ba.size() != flev + 1) {
      ba.resize(flev + 1);
      dm.resize(flev + 1);
      changed = true;
    } else {
      for (int i = 0; i <= flev && !changed; i++) {
        // Check if BoxArrays have changed during regridding
        if (ba[i] != this->boxArray(i)) {
          changed = true;
        }

        if (!changed) {
          // Check DistributionMaps have changed during regridding
          if (dm[i] != this->DistributionMap(i)) {
            changed = true;
          }
        }
      }
    }

    if (changed) {
      //
      // We only need to call Redistribute if the BoxArrays or DistMaps have
      // changed.
      //
      if (spray_verbose && ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Calling redistribute because grid has changed "
                       << '\n';
      }
      if (flev == 0) {
        // Do a local redistribute
        theSprayPC()->Redistribute(lbase, theSprayPC()->finestLevel(), 1);
      } else {
        theSprayPC()->Redistribute(lbase, theSprayPC()->finestLevel(), 1);
      }
      //
      // Use the new BoxArray and DistMap to define ba and dm for next time.
      //
      for (int i = 0; i <= flev; i++) {
        ba[i] = this->boxArray(i);
        dm[i] = this->DistributionMap(i);
      }
    } else {
      if (verbose && ParallelDescriptor::IOProcessor()) {
        amrex::Print()
          << "NOT calling redistribute because grid has NOT changed " << '\n';
      }
    }
  }
}

#endif
