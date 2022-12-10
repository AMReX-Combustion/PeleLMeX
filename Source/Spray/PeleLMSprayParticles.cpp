
#include <PeleLM.H>

#ifdef PELELM_USE_SPRAY
#include "SprayParticles.H"
#include <AMReX_FillPatchUtil.H>

using namespace amrex;

namespace {
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
Vector<Real> spray_cfl;
Vector<int> spray_state_ghosts;
Vector<int> spray_source_ghosts;
Vector<int> spray_ghost_num;
Vector<int> prev_state;
Vector<int> prev_source;

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
bool mesh_regrid = true;
std::string init_file;
int init_function = 1;
int spray_verbose = 0;
Real max_spray_cfl = 5.;
} // namespace

bool PeleLM::do_spray_particles = true;
// momentum + density + fuel species + enthalpy
int PeleLM::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;

int PeleLM::write_spray_ascii_files = 0;
int PeleLM::plot_spray_src = 0;

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
PeleLM::SprayEstDt()
{
  Real estdt = 1.0e200;
  if (!do_spray_particles || theSprayPC() == nullptr) {
    return estdt;
  }
  BL_PROFILE("PeleLM::SprayEstDt()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    Real estdt_lev = theSprayPC()->estTimestep(lev, max_spray_cfl);
    if (estdt_lev > 0. && estdt_lev < estdt) {
      estdt = estdt_lev;
    } else if (lev > 0) {
      // Limit time step as if particles are on finest level
      estdt /= Real(refRatio(lev - 1)[0]);
    }
  }
  return estdt;
}

void
PeleLM::SprayReadParameters()
{
  ParmParse pp("peleLM");

  pp.query("do_spray_particles", do_spray_particles);
  if (!do_spray_particles) {
    return;
  }
  const Real temp_cfl = max_spray_cfl;
  // Mush change dtmod to 1 since we only do MKD
  sprayData.dtmod = 1.;
  SprayParticleContainer::readSprayParams(
    spray_verbose, max_spray_cfl, write_spray_ascii_files,
    plot_spray_src, init_function, init_file, sprayData, temp_cfl);
}

void
PeleLM::SpraySetup()
{
  if (!do_spray_particles) {
    return;
  }
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    amrex::Abort("Cannot have more spray fuel species than fluid species");
  }
  SprayParticleContainer::spraySetup(sprayData);
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
  if (theSprayPC() != nullptr) {
    if (level < finest_level) {
      SprayParticleContainer::AoS virts;
      setupVirtualParticles(level + 1);
      theVirtPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);

      theSprayPC()->CreateVirtualParticles(level + 1, virts);
      theVirtPC()->AddParticlesAtLevel(virts, level);
    }
  }
}

void
PeleLM::removeVirtualParticles(const int level)
{
  if (theVirtPC() != nullptr) {
    theVirtPC()->RemoveParticlesAtLevel(level);
  }
}

void
PeleLM::setupGhostParticles(const int ngrow, const int level)
{
  BL_PROFILE("PeleLM::setupGhostParticles()");
  AMREX_ASSERT(level < finest_level);
  if (theSprayPC() != nullptr) {
    SprayParticleContainer::AoS ghosts;
    theSprayPC()->CreateGhostParticles(level, ngrow, ghosts);
    theGhostPC()->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleLM::removeGhostParticles(const int level)
{
  if (GhostPC != nullptr) {
    GhostPC->RemoveParticlesAtLevel(level);
  }
}

void
PeleLM::SprayCreateData()
{
  SprayPC = new SprayParticleContainer(
    this, &m_phys_bc, sprayData, scomps, max_spray_cfl);
  theSprayPC()->SetVerbose(spray_verbose);
  VirtPC = new SprayParticleContainer(
    this, &m_phys_bc, sprayData, scomps, max_spray_cfl);
  GhostPC = new SprayParticleContainer(
    this, &m_phys_bc, sprayData, scomps, max_spray_cfl);
}

void
PeleLM::SprayInit()
{
  BL_PROFILE("PeleLM::SprayInit()");

  //
  // Make sure to call RemoveParticlesOnExit() on exit.
  //
  amrex::ExecOnFinalize(RemoveParticlesOnExit);

  if (!do_spray_particles) {
    return;
  }
  SprayCreateData();
  // Check if we are restarting and if we restart file has particles
  std::string restart_file =
    (m_restart_chkfile.empty()) ? m_restart_pltfile : m_restart_chkfile;
  std::string restart_partfile = restart_file + "/particles";
  if (FileSystem::Exists(restart_partfile)) {
    SprayRestart(restart_file);
  } else {
    if (!m_restart_chkfile.empty() || !m_restart_pltfile.empty()) {
      std::string warn_msg = "Restart file does not contain particles. "
                             "Particles are being initialized from scratch.";
      amrex::Warning(warn_msg);
    }

    bool init_part = true;
    if (!init_file.empty()) {
      theSprayPC()->InitFromAsciiFile(init_file, NSR_SPR + NAR_SPR);
    }
    if (init_function <= 0) {
      init_part = false;
    }
    ProbParm const* lprobparm = prob_parm;
    theSprayPC()->InitSprayParticles(init_part, *lprobparm);
    SprayPostRegrid();
    SprayInjectRedist();
    if (spray_verbose >= 1) {
      amrex::Print() << "Total number of initial particles "
                     << theSprayPC()->TotalNumberOfParticles(false, false)
                     << std::endl;
    }
  }
}

void
PeleLM::SprayRestart(const std::string& restart_file)
{
  if (do_spray_particles) {
    //
    // Make sure to call RemoveParticlesOnExit() on exit.
    //
    amrex::ExecOnFinalize(RemoveParticlesOnExit);
    {
      theSprayPC()->Restart(restart_file, "particles");
      ProbParm const* lprobparm = prob_parm;
      theSprayPC()->InitSprayParticles(false, *lprobparm);
      SprayPostRegrid();
      amrex::Gpu::Device::streamSynchronize();
    }
  }
}

// Sets the number of ghost cells for the state, source, and ghost particles
// Also creates and fills the state data
void
PeleLM::SpraySetState(const Real& a_flow_dt)
{
  if (!do_spray_particles) {
    return;
  }
  spray_cfl.resize(finest_level + 1);
  // Number of ghost cells needed to interpolate state to spray position
  spray_state_ghosts.resize(finest_level + 1);
  // Number of ghost cells needed to interpolate spray source to mesh
  spray_source_ghosts.resize(finest_level + 1);
  // Number of ghost cells in which to make ghost particles
  spray_ghost_num.resize(finest_level + 1);
  // Determine the max velocity of particles at each level
  // This is necessary because not every level will be aware of the fastest
  // particle for when ghost and virtual particles are present
  Real max_vel = 0.;
  spray_state_ghosts[0] = 0;
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto const dx = geom[lev].CellSizeArray();
    // Extract velocity and CFL from a given spray CFL
    Real spraydt_lev = theSprayPC()->estTimestep(lev, max_spray_cfl);
    Real vel_lev = max_spray_cfl * dx[0] / spraydt_lev;
    max_vel = amrex::max(max_vel, vel_lev);
    if (spraydt_lev > 0.) {
      spray_cfl[lev] = max_spray_cfl / spraydt_lev * a_flow_dt;
    } else {
      spray_cfl[lev] = 0.;
    }
  }
  for (int lev = finest_level; lev >= 0; --lev) {
    auto const dx = geom[lev].CellSizeArray();
    spray_cfl[lev] = amrex::max(spray_cfl[lev], max_vel * a_flow_dt / dx[0]);
    if (lev < finest_level) {
      // Note: Ghost particles made during level N depend on information at
      // level N+1
      spray_ghost_num[lev] = SprayParticleContainer::getGhostPartCells(
        lev + 1, finest_level, 1, spray_cfl[lev + 1]);
    } else {
      spray_ghost_num[lev] = 0;
    }
    int state_ghosts = SprayParticleContainer::getStateGhostCells(
      lev, finest_level, 1, spray_cfl[lev]);
    int source_ghosts = SprayParticleContainer::getSourceGhostCells(
      lev, finest_level, 1, spray_cfl[lev]);
    spray_state_ghosts[lev] = state_ghosts;
    spray_source_ghosts[lev] = source_ghosts;
    if (
      mesh_regrid || prev_state[lev] != state_ghosts ||
      prev_source[lev] != source_ghosts) {
      m_spraystate[lev].reset(new MultiFab(
        grids[lev], dmap[lev], NVAR, state_ghosts, MFInfo(), *m_factory[lev]));
      m_spraysource[lev].reset(new MultiFab(
        grids[lev], dmap[lev], num_spray_src, source_ghosts, MFInfo(),
        *m_factory[lev]));
    }
    fillpatch_state(lev, m_cur_time, *(m_spraystate[lev].get()), state_ghosts);
    m_spraysource[lev]->setVal(0.);
  }
  mesh_regrid = false;
}

void
PeleLM::SprayAddSource(const int level)
{
  MultiFab& source = *(m_spraysource[level].get());
  MultiFab& extsource = *(m_extSource[level].get());
  const int eghosts = extsource.nGrow();
  MultiFab::Add(
    extsource, source, scomps.rhoSrcIndx, scomps.rhoIndx, 1, eghosts);
  MultiFab::Add(
    extsource, source, scomps.engSrcIndx, scomps.engIndx, 1, eghosts);
  MultiFab::Add(
    extsource, source, scomps.momSrcIndx, scomps.momIndx, AMREX_SPACEDIM,
    eghosts);
  for (int n = 0; n < SPRAY_FUEL_NUM; ++n) {
    const int dstcomp = scomps.specIndx + sprayData.dep_indx[n];
    MultiFab::Add(
      extsource, source, scomps.specSrcIndx + n, dstcomp, 1, eghosts);
  }
}

void
PeleLM::SprayMKD(const Real time, const Real dt)
{
  if (!do_spray_particles) {
    return;
  }
  if (spray_verbose) {
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";
  }
  BL_PROFILE("PeleLM::SprayMKD()");
  // Setup the virtual particles that represent particles on finer levels
  setupVirtualParticles(0);
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (spray_verbose > 1) {
      amrex::Print() << "SprayMKDLevel " << lev << std::endl;
    }
    SprayMKDLevel(lev, time, dt);
    SprayAddSource(lev);
    removeGhostParticles(lev);
    removeVirtualParticles(lev);
    m_spraysource[lev]->setVal(0.);
  }
}

void
PeleLM::SprayMKDLevel(const int level, const Real time, const Real dt)
{
  if (level < finest_level) {
    // Make a copy of the particles on this level into ghost particles
    // for the finer level
    int ghost_width = spray_ghost_num[level];
    setupGhostParticles(ghost_width, level);
  }

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  auto const* ltransparm = PeleLM::trans_parms.device_trans_parm();

  amrex::MultiFab& state = *(m_spraystate[level].get());
  amrex::MultiFab& source = *(m_spraysource[level].get());
  const int state_ghosts = spray_state_ghosts[level];
  const int source_ghosts = spray_source_ghosts[level];
  bool isVirt = false;
  bool isGhost = false;
  bool doMove = true;
  // Do the valid particles themselves
  theSprayPC()->moveKickDrift(
    state, source, level, dt, time, isVirt, isGhost, state_ghosts,
    source_ghosts, doMove, ltransparm, spray_cfl[level]);
  if (level < finest_level) {
    isVirt = true;
    isGhost = false;
    theVirtPC()->moveKickDrift(
      state, source, level, dt, time, isVirt, isGhost, state_ghosts,
      source_ghosts, doMove, ltransparm, spray_cfl[level]);
  }
  if (theGhostPC() != nullptr && level != 0) {
    isVirt = false;
    isGhost = true;
    theGhostPC()->moveKickDrift(
      state, source, level, dt, time, isVirt, isGhost, state_ghosts,
      source_ghosts, doMove, ltransparm, spray_cfl[level]);
  }
  source.SumBoundary(geom[level].periodicity());
}

void
PeleLM::SprayPostRegrid()
{
  static Vector<BoxArray> ba_spray;
  static Vector<DistributionMapping> dm_spray;
  bool changed = false;
  if (ba_spray.size() != finest_level + 1) {
    ba_spray.resize(finest_level + 1);
    dm_spray.resize(finest_level + 1);
    prev_state.resize(finest_level + 1);
    prev_source.resize(finest_level + 1);
    changed = true;
  } else {
    for (int lev = 0; lev <= finest_level && !changed; lev++) {
      if (ba_spray[lev] != grids[lev]) {
        changed = true;
      }
      if (!changed && dm_spray[lev] != dmap[lev]) {
        changed = true;
      }
    }
  }
  // Update the local BoxArray and DistributionMap
  if (changed) {
    mesh_regrid = true;
    for (int lev = 0; lev <= finest_level; ++lev) {
      ba_spray[lev] = grids[lev];
      dm_spray[lev] = dmap[lev];
      prev_state[lev] = -1;
      prev_source[lev] = -1;
    }
    theSprayPC()->Redistribute();
  }
}

void
PeleLM::SprayInjectRedist()
{
  BL_PROFILE("PeleLM::SprayInjectRedist");
  Long prev_count = 0;
  if (spray_verbose >= 3) {
    prev_count = theSprayPC()->TotalNumberOfParticles(true, false);
  }
  bool injected = false;
  for (int lev = 0; lev <= finest_level; ++lev) {
    int nstep = 0; // Unused
    ProbParm const* lprobparm = prob_parm;
    Real cur_time = m_t_new[lev]; // Still the time from the last time step
    Real dt = m_dt;
    bool lev_injected = theSprayPC()->injectParticles(
      cur_time, dt, nstep, lev, finest_level, *lprobparm);
    if (lev_injected) {
      injected = true;
    }
  }
  // We must redistribute after each time step
  theSprayPC()->Redistribute();
  if (spray_verbose >= 3 && injected) {
    Long new_count = theSprayPC()->TotalNumberOfParticles(true, false);
    Long num_inj = new_count - prev_count;
    amrex::Print() << "Injected " << num_inj << " particles at time "
                   << m_t_new[0] << std::endl;
  }
}

#endif
