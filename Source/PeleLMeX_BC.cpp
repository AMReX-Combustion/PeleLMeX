#include <PeleLMeX.H>
#include <PeleLMeX_BCfill.H>
#include <AMReX_FillPatchUtil.H>
#include <memory>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#include <AMReX_EB2.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFabFactory.H>
#endif

// Conversion from physBC to fieldBC maps
// Components are  Interior, Inflow, Outflow, Symmetry, &
// SlipWallAdiab, NoSlipWallAdiab, SlipWallIsoTherm, NoSlipWallIsoTherm.

constexpr int norm_vel_bc[] = {
  amrex::BCType::int_dir,     amrex::BCType::ext_dir, amrex::BCType::foextrap,
  amrex::BCType::reflect_odd, amrex::BCType::ext_dir, amrex::BCType::ext_dir,
  amrex::BCType::ext_dir,     amrex::BCType::ext_dir};

constexpr int tang_vel_bc[] = {
  amrex::BCType::int_dir,      amrex::BCType::ext_dir,  amrex::BCType::foextrap,
  amrex::BCType::reflect_even, amrex::BCType::hoextrap, amrex::BCType::ext_dir,
  amrex::BCType::hoextrap,     amrex::BCType::ext_dir};

constexpr int density_bc[] = {
  amrex::BCType::int_dir,      amrex::BCType::ext_dir,  amrex::BCType::foextrap,
  amrex::BCType::reflect_even, amrex::BCType::foextrap, amrex::BCType::foextrap,
  amrex::BCType::foextrap,     amrex::BCType::foextrap};

constexpr int species_bc[] = {
  amrex::BCType::int_dir,      amrex::BCType::ext_dir,  amrex::BCType::foextrap,
  amrex::BCType::reflect_even, amrex::BCType::foextrap, amrex::BCType::foextrap,
  amrex::BCType::foextrap,     amrex::BCType::foextrap};

constexpr int ode_bc[] = {amrex::BCType::int_dir,  amrex::BCType::ext_dir,
                          amrex::BCType::foextrap, amrex::BCType::reflect_even,
                          amrex::BCType::foextrap, amrex::BCType::foextrap,
                          amrex::BCType::foextrap, amrex::BCType::foextrap};

constexpr int rhoh_bc[] = {amrex::BCType::int_dir,  amrex::BCType::ext_dir,
                           amrex::BCType::foextrap, amrex::BCType::reflect_even,
                           amrex::BCType::foextrap, amrex::BCType::foextrap,
                           amrex::BCType::ext_dir,  amrex::BCType::ext_dir};

constexpr int temp_bc[] = {amrex::BCType::int_dir,  amrex::BCType::ext_dir,
                           amrex::BCType::foextrap, amrex::BCType::reflect_even,
                           amrex::BCType::foextrap, amrex::BCType::foextrap,
                           amrex::BCType::ext_dir,  amrex::BCType::ext_dir};

constexpr int aux_bc[] = {amrex::BCType::int_dir,  amrex::BCType::ext_dir,
                          amrex::BCType::foextrap, amrex::BCType::reflect_even,
                          amrex::BCType::foextrap, amrex::BCType::foextrap,
                          amrex::BCType::ext_dir,  amrex::BCType::ext_dir};

constexpr int divu_bc[] = {
  amrex::BCType::int_dir,      amrex::BCType::reflect_even,
  amrex::BCType::reflect_even, amrex::BCType::reflect_even,
  amrex::BCType::reflect_even, amrex::BCType::reflect_even,
  amrex::BCType::reflect_even, amrex::BCType::reflect_even};

// Following incflo rather than IAMR here
constexpr int force_bc[] = {amrex::BCType::int_dir,  amrex::BCType::foextrap,
                            amrex::BCType::foextrap, amrex::BCType::foextrap,
                            amrex::BCType::foextrap, amrex::BCType::foextrap,
                            amrex::BCType::foextrap, amrex::BCType::foextrap};

#ifdef PELE_USE_PLASMA
constexpr int nE_bc[] = {
  amrex::BCType::int_dir,      amrex::BCType::ext_dir,
  amrex::BCType::foextrap,     amrex::BCType::reflect_even,
  amrex::BCType::reflect_even, amrex::BCType::reflect_even,
  amrex::BCType::ext_dir,      amrex::BCType::ext_dir};

constexpr int phiV_bc[] = {
  amrex::BCType::int_dir, amrex::BCType::ext_dir, amrex::BCType::reflect_even};
#endif

#ifdef PELE_USE_SOOT
constexpr int soot_bc[] = {
  amrex::BCType::int_dir,      amrex::BCType::ext_dir,
  amrex::BCType::foextrap,     amrex::BCType::reflect_even,
  amrex::BCType::reflect_even, amrex::BCType::reflect_even,
  amrex::BCType::ext_dir,      amrex::BCType::ext_dir};
#endif

amrex::InterpBase*
PeleLM::
  getInterpolator( // NOLINT(readability-convert-member-functions-to-static)
    const int a_method) const
{
  amrex::InterpBase* mapper = nullptr;

  switch (a_method) {
  case 0:
    mapper = &amrex::mf_pc_interp;
    break;

  case 1:
#ifdef AMREX_USE_EB
    // Get EB-aware interpolator when needed
    mapper = (EBFactory(0).isAllRegular()) ? &amrex::mf_cell_cons_interp
                                           : &amrex::eb_mf_cell_cons_interp;
#else
    mapper = &amrex::mf_cell_cons_interp;
#endif
    break;

  case 2:
#ifdef AMREX_USE_EB
    amrex::Abort("Regrid interpolation method = 2 not available with EB !");
#else
    mapper = &amrex::mf_linear_slope_minmax_interp;
#endif
    break;

  default:
    amrex::Abort("Unknown interpolation method");
  }

  return mapper;
}

void
PeleLM::setBoundaryConditions()
{

  // Initialize the BCRecs
  m_bcrec_state.resize(NVAR);
  constexpr int sizeForceBC = amrex::max<int>(AMREX_SPACEDIM, NUM_SPECIES + 2);
  m_bcrec_force.resize(sizeForceBC);
  m_bcrec_aux.resize(m_nAux);

  // Convert m_phys_bc into field BCs
  // Get m_phys_bc
  const int* lo_bc = m_phys_bc.lo();
  const int* hi_bc = m_phys_bc.hi();

  // Velocity
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    for (int idim2 = 0; idim2 < AMREX_SPACEDIM; ++idim2) {
      if (idim == idim2) {
        m_bcrec_state[VELX + idim].setLo(idim2, norm_vel_bc[lo_bc[idim2]]);
        m_bcrec_state[VELX + idim].setHi(idim2, norm_vel_bc[hi_bc[idim2]]);
      } else {
        m_bcrec_state[VELX + idim].setLo(idim2, tang_vel_bc[lo_bc[idim2]]);
        m_bcrec_state[VELX + idim].setHi(idim2, tang_vel_bc[hi_bc[idim2]]);
      }
    }
  }

  // General forces: use int_dir in interior and foextrap otherwise
  for (int i = 0; i < sizeForceBC; ++i) {
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_force[i].setLo(idim, force_bc[lo_bc[idim]]);
      m_bcrec_force[i].setHi(idim, force_bc[hi_bc[idim]]);
    }
  }

  if (m_incompressible == 0) {
    // Density
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_state[DENSITY].setLo(idim, density_bc[lo_bc[idim]]);
      m_bcrec_state[DENSITY].setHi(idim, density_bc[hi_bc[idim]]);
    }

    // Species
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        m_bcrec_state[FIRSTSPEC + n].setLo(idim, species_bc[lo_bc[idim]]);
        m_bcrec_state[FIRSTSPEC + n].setHi(idim, species_bc[hi_bc[idim]]);
      }
    }

    // Enthalpy
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_state[RHOH].setLo(idim, rhoh_bc[lo_bc[idim]]);
      m_bcrec_state[RHOH].setHi(idim, rhoh_bc[hi_bc[idim]]);
    }

    // Temperature
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_state[TEMP].setLo(idim, temp_bc[lo_bc[idim]]);
      m_bcrec_state[TEMP].setHi(idim, temp_bc[hi_bc[idim]]);
    }

    // rhoRT: reflect even on all but interior bndy
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_state[RHORT].setLo(idim, divu_bc[lo_bc[idim]]);
      m_bcrec_state[RHORT].setHi(idim, divu_bc[hi_bc[idim]]);
    }

    // ODEs
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      for (int n = 0; n < NUM_ODE; ++n) {
        m_bcrec_state[FIRSTODE + n].setLo(idim, ode_bc[lo_bc[idim]]);
        m_bcrec_state[FIRSTODE + n].setHi(idim, ode_bc[hi_bc[idim]]);
      }
    }

    // divU
    if (m_has_divu != 0) {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_bcrec_divu.setLo(idim, divu_bc[lo_bc[idim]]);
        m_bcrec_divu.setHi(idim, divu_bc[hi_bc[idim]]);
      }
    }
    // auxiliaries - assumed to be the same as species
    for (int n = 0; n < m_nAux; ++n) {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        m_bcrec_aux[n].setLo(idim, aux_bc[lo_bc[idim]]);
        m_bcrec_aux[n].setHi(idim, aux_bc[hi_bc[idim]]);
      }
    }

#ifdef PELE_USE_PLASMA
    // nE
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_state[NE].setLo(idim, nE_bc[lo_bc[idim]]);
      m_bcrec_state[NE].setHi(idim, nE_bc[hi_bc[idim]]);
    }

    // Get m_phiV_bc
    const int* lo_phibc = m_phiV_bc.lo();
    const int* hi_phibc = m_phiV_bc.hi();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_state[PHIV].setLo(idim, phiV_bc[lo_phibc[idim]]);
      m_bcrec_state[PHIV].setHi(idim, phiV_bc[hi_phibc[idim]]);
    }

    // Hack charged species BCs
    int FIRSTIONinVar = FIRSTSPEC + NUM_SPECIES - NUM_IONS;
    int FIRSTIONinSpec = NUM_SPECIES - NUM_IONS;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      for (int n = 0; n < NUM_IONS; ++n) {
        auto const bcIonSave = m_bcrec_state[FIRSTIONinVar + n];
        m_bcrec_state[FIRSTIONinVar + n] =
          hackBCChargedParticle(zk[FIRSTIONinSpec + n], bcIonSave);
      }
    }
    // Need to hack nE too actually ...
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      auto const bcnESave = m_bcrec_state[NE];
      m_bcrec_state[NE] = hackBCChargedParticle(-1.0, bcnESave);
    }
#endif
#ifdef PELE_USE_SOOT
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      for (int mom = 0; mom < NUMSOOTVAR; ++mom) {
        m_bcrec_state[FIRSTSOOT + mom].setLo(idim, soot_bc[lo_bc[idim]]);
        m_bcrec_state[FIRSTSOOT + mom].setHi(idim, soot_bc[hi_bc[idim]]);
      }
    }
#endif
    // Dummy
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      m_bcrec_dummy.setLo(idim, amrex::BCType::int_dir);
      m_bcrec_dummy.setHi(idim, amrex::BCType::int_dir);
    }
  }
}

amrex::Vector<amrex::BCRec>
PeleLM::fetchBCRecArray(const int scomp, const int ncomp)
{
  amrex::Vector<amrex::BCRec> bc(ncomp);
  for (int comp = 0; comp < ncomp; ++comp) {
    bc[comp] = m_bcrec_state[scomp + comp];
  }
  return bc;
}

amrex::Vector<amrex::BCRec>
PeleLM::fetchBCRecAuxArray(const int scomp, const int ncomp)
{
  amrex::Vector<amrex::BCRec> bc(ncomp);
  for (int comp = 0; comp < ncomp; ++comp) {
    bc[comp] = m_bcrec_aux[scomp + comp];
  }
  return bc;
}

amrex::Vector<amrex::BCRec>
PeleLM::fetchBCRecDummyArray(const int /* scomp */, const int ncomp)
{
  amrex::Vector<amrex::BCRec> bc(ncomp);
  for (int comp = 0; comp < ncomp; ++comp) {
    bc[comp] = m_bcrec_dummy;
  }
  return bc;
}

//-----------------------------------------------------------------------------
// The following work directly on the leveldata
// Fill the entire class state at once
void
PeleLM::fillPatchState(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchState()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    fillPatchState(lev, a_time);
  }
}

// Fill the a given level class state
void
PeleLM::fillPatchState(const int lev, const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchStateLev()");

  auto* ldata_p = getLevelDataPtr(lev, a_time);
  const amrex::Real time = getTime(lev, a_time);

  fillpatch_state(lev, time, ldata_p->state, m_nGrowState);
  if (m_incompressible == 0) {
    if (m_has_divu != 0) {
      fillpatch_divu(lev, time, ldata_p->divu, ldata_p->divu.nGrow());
    }
  }
}

// Fill a state components
void
PeleLM::fillPatchDensity(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchDensity()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto* ldata_p = getLevelDataPtr(lev, a_time);
    const amrex::Real time = getTime(lev, a_time);
    fillpatch_density(lev, time, ldata_p->state, DENSITY, m_nGrowState);
  }
}

void
PeleLM::fillPatchSpecies(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchSpecies()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto* ldata_p = getLevelDataPtr(lev, a_time);
    const amrex::Real time = getTime(lev, a_time);
    fillpatch_species(lev, time, ldata_p->state, FIRSTSPEC, m_nGrowState);
  }
}

void
PeleLM::fillPatchTemp(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchTemp()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto* ldata_p = getLevelDataPtr(lev, a_time);
    const amrex::Real time = getTime(lev, a_time);
    fillpatch_temp(lev, time, ldata_p->state, TEMP, m_nGrowState);
  }
}

void
PeleLM::fillPatchAux(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchAux()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto* ldata_p = getLevelDataPtr(lev, a_time);
    const amrex::Real time = getTime(lev, a_time);
    fillpatch_aux(lev, time, ldata_p->auxiliaries, m_nGrowState);
  }
}

#ifdef PELE_USE_PLASMA
void
PeleLM::fillPatchPhiV(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::fillPatchPhiV()");
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto ldata_p = getLevelDataPtr(lev, a_time);
    const amrex::Real time = getTime(lev, a_time);
    fillpatch_phiV(lev, time, ldata_p->state, PHIV, m_nGrowState);
  }
}
#endif
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// The following return a fillpatched MF ptr at a given level
// Fill the entire state at once
std::unique_ptr<amrex::MultiFab>
PeleLM::fillPatchState(const int lev, const amrex::Real a_time, const int nGrow)
{
  BL_PROFILE("PeleLMeX::fillPatchState()");

  std::unique_ptr<amrex::MultiFab> mf;
  if (m_incompressible != 0) {
    mf = std::make_unique<amrex::MultiFab>(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow, amrex::MFInfo(),
      Factory(lev));
  } else {
    mf = std::make_unique<amrex::MultiFab>(
      grids[lev], dmap[lev], NVAR, nGrow, amrex::MFInfo(), Factory(lev));
  }
  fillpatch_state(lev, a_time, *mf, nGrow);

  return mf;
}

std::unique_ptr<amrex::MultiFab>
PeleLM::fillPatchReact(const int lev, const amrex::Real a_time, const int nGrow)
{
  BL_PROFILE("PeleLMeX::fillPatchReact()");

#ifdef PELE_USE_PLASMA
  constexpr int IRsize = NUM_SPECIES + 1;
#else
  constexpr int IRsize = NUM_SPECIES;
#endif
  std::unique_ptr<amrex::MultiFab> mf;
  mf = std::make_unique<amrex::MultiFab>(
    grids[lev], dmap[lev], IRsize, nGrow, amrex::MFInfo(), Factory(lev));
  fillpatch_reaction(lev, a_time, *mf, nGrow);

  return mf;
}
//-----------------------------------------------------------------------------

// Fill the state
void
PeleLM::fillpatch_state(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_state,
  const int nGhost)
{
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();

  const int nCompState = (m_incompressible) != 0 ? AMREX_SPACEDIM : NVAR;

  fillTurbInflow(a_state, VELX, lev, a_time);

  if (m_use_inlet_from_plane != 0) {
    fillFromRecyclingPlane(a_state, 0, lev);
  }

  if (lev == 0) {
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
      bndry_func(
        geom[lev], fetchBCRecArray(0, nCompState),
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux,
          static_cast<int>(turb_inflow.is_initialized()),
          static_cast<int>(m_use_inlet_from_plane)});
    FillPatchSingleLevel(
      a_state, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, nCompState, geom[lev], bndry_func, 0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
      crse_bndry_func(
        geom[lev - 1], fetchBCRecArray(0, nCompState),
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux,
          static_cast<int>(turb_inflow.is_initialized()),
          static_cast<int>(m_use_inlet_from_plane)});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
      fine_bndry_func(
        geom[lev], fetchBCRecArray(0, nCompState),
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux,
          static_cast<int>(turb_inflow.is_initialized()),
          static_cast<int>(m_use_inlet_from_plane)});
    FillPatchTwoLevels(
      a_state, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->state), &(m_leveldata_new[lev - 1]->state)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, nCompState, geom[lev - 1], geom[lev],
      crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1), mapper,
      fetchBCRecArray(0, nCompState), 0);
  }

  a_state.EnforcePeriodicity(geom[lev].periodicity());
}

// Fill the density
void
PeleLM::fillpatch_density(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_density,
  const int rho_comp,
  const int nGhost)
{
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  if (lev == 0) {

    // Density
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDens<ProblemSpecificFunctions>>>
      bndry_func_rho(
        geom[lev], fetchBCRecArray(DENSITY, 1),
        PeleLMCCFillExtDirDens<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchSingleLevel(
      a_density, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, DENSITY, rho_comp, 1, geom[lev],
      bndry_func_rho, 0);

  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    // Density
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDens<ProblemSpecificFunctions>>>
      crse_bndry_func_rho(
        geom[lev - 1], fetchBCRecArray(DENSITY, 1),
        PeleLMCCFillExtDirDens<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDens<ProblemSpecificFunctions>>>
      fine_bndry_func_rho(
        geom[lev], fetchBCRecArray(DENSITY, 1),
        PeleLMCCFillExtDirDens<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchTwoLevels(
      a_density, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->state), &(m_leveldata_new[lev - 1]->state)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, DENSITY, rho_comp, 1, geom[lev - 1],
      geom[lev], crse_bndry_func_rho, 0, fine_bndry_func_rho, 0,
      refRatio(lev - 1), mapper, fetchBCRecArray(DENSITY, 1), 0);
  }
}

// Fill the mass fractions
void
PeleLM::fillpatch_species(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_species,
  const int rhoY_comp,
  const int nGhost)
{
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  if (lev == 0) {

    // Species
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirSpec<ProblemSpecificFunctions>>>
      bndry_func(
        geom[lev], fetchBCRecArray(FIRSTSPEC, NUM_SPECIES),
        PeleLMCCFillExtDirSpec<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchSingleLevel(
      a_species, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, FIRSTSPEC, rhoY_comp, NUM_SPECIES,
      geom[lev], bndry_func, 0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    // Species
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirSpec<ProblemSpecificFunctions>>>
      crse_bndry_func(
        geom[lev - 1], fetchBCRecArray(FIRSTSPEC, NUM_SPECIES),
        PeleLMCCFillExtDirSpec<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirSpec<ProblemSpecificFunctions>>>
      fine_bndry_func(
        geom[lev], fetchBCRecArray(FIRSTSPEC, NUM_SPECIES),
        PeleLMCCFillExtDirSpec<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchTwoLevels(
      a_species, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->state), &(m_leveldata_new[lev - 1]->state)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, FIRSTSPEC, rhoY_comp, NUM_SPECIES,
      geom[lev - 1], geom[lev], crse_bndry_func, 0, fine_bndry_func, 0,
      refRatio(lev - 1), mapper, fetchBCRecArray(FIRSTSPEC, NUM_SPECIES), 0);
  }
}

// Fill temperature
void
PeleLM::fillpatch_temp(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_temp,
  const int temp_comp,
  const int nGhost)
{
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  if (lev == 0) {
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirTemp<ProblemSpecificFunctions>>>
      bndry_func(
        geom[lev], fetchBCRecArray(TEMP, 1),
        PeleLMCCFillExtDirTemp<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchSingleLevel(
      a_temp, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, TEMP, temp_comp, 1, geom[lev], bndry_func,
      0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirTemp<ProblemSpecificFunctions>>>
      crse_bndry_func(
        geom[lev - 1], fetchBCRecArray(TEMP, 1),
        PeleLMCCFillExtDirTemp<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirTemp<ProblemSpecificFunctions>>>
      fine_bndry_func(
        geom[lev], fetchBCRecArray(TEMP, 1),
        PeleLMCCFillExtDirTemp<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchTwoLevels(
      a_temp, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->state), &(m_leveldata_new[lev - 1]->state)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, TEMP, temp_comp, 1, geom[lev - 1],
      geom[lev], crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1),
      mapper, fetchBCRecArray(TEMP, 1), 0);
  }
}

// Fill the auxiliaries
void
PeleLM::fillpatch_aux(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_aux,
  const int nGhost)
{

  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();

  if (lev == 0) {
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirAux<ProblemSpecificFunctions>>>
      bndry_func(
        geom[lev], fetchBCRecAuxArray(0, m_nAux),
        PeleLMCCFillExtDirAux<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchSingleLevel(
      a_aux, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->auxiliaries),
       &(m_leveldata_new[lev]->auxiliaries)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, m_nAux, geom[lev], bndry_func, 0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirAux<ProblemSpecificFunctions>>>
      crse_bndry_func(
        geom[lev - 1], fetchBCRecAuxArray(0, m_nAux),
        PeleLMCCFillExtDirAux<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirAux<ProblemSpecificFunctions>>>
      fine_bndry_func(
        geom[lev], fetchBCRecAuxArray(0, m_nAux),
        PeleLMCCFillExtDirAux<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchTwoLevels(
      a_aux, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->auxiliaries),
       &(m_leveldata_new[lev - 1]->auxiliaries)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->auxiliaries),
       &(m_leveldata_new[lev]->auxiliaries)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, m_nAux, geom[lev - 1], geom[lev],
      crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1), mapper,
      fetchBCRecAuxArray(0, m_nAux), 0);
  }

  a_aux.EnforcePeriodicity(geom[lev].periodicity());
}

#ifdef PELE_USE_PLASMA
// Fill electro-static potential
void
PeleLM::fillpatch_phiV(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_temp,
  const int phiV_comp,
  const int nGhost)
{
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  if (lev == 0) {
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirPhiV<ProblemSpecificFunctions>>>
      bndry_func(
        geom[lev], fetchBCRecArray(PHIV, 1),
        PeleLMCCFillExtDirPhiV<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchSingleLevel(
      a_temp, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, PHIV, phiV_comp, 1, geom[lev], bndry_func,
      0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirPhiV<ProblemSpecificFunctions>>>
      crse_bndry_func(
        geom[lev - 1], fetchBCRecArray(PHIV, 1),
        PeleLMCCFillExtDirPhiV<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirPhiV<ProblemSpecificFunctions>>>
      fine_bndry_func(
        geom[lev], fetchBCRecArray(PHIV, 1),
        PeleLMCCFillExtDirPhiV<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux});
    FillPatchTwoLevels(
      a_temp, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->state), &(m_leveldata_new[lev - 1]->state)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->state), &(m_leveldata_new[lev]->state)},
      {m_t_old[lev], m_t_new[lev]}, PHIV, phiV_comp, 1, geom[lev - 1],
      geom[lev], crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1),
      mapper, fetchBCRecArray(PHIV, 1), 0);
  }
}
#endif

// Fill the divU
void
PeleLM::fillpatch_divu(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_divu,
  const int nGhost)
{
  if (lev == 0) {
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      bndry_func(geom[lev], {m_bcrec_divu}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchSingleLevel(
      a_divu, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->divu), &(m_leveldata_new[lev]->divu)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev], bndry_func, 0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      crse_bndry_func(
        geom[lev - 1], {m_bcrec_divu}, PeleLMCCFillExtDirDummy{m_nAux});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      fine_bndry_func(
        geom[lev], {m_bcrec_divu}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchTwoLevels(
      a_divu, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->divu), &(m_leveldata_new[lev - 1]->divu)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->divu), &(m_leveldata_new[lev]->divu)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, 1, geom[lev - 1], geom[lev],
      crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1), mapper,
      {m_bcrec_divu}, 0);
  }
}

// Fillpatch a vector of forces:
// -> actually only modifies the ghost cells : fillBoundary, C/F interp,
// foextrap on domain BCs
void
PeleLM::fillpatch_forces(
  const amrex::Real a_time,
  amrex::Vector<amrex::MultiFab*> const& a_force,
  const int nGrowForce)
{
  AMREX_ASSERT(a_force[0]->nComp() <= m_bcrec_force.size());
  const int nComp = a_force[0]->nComp();

  int lev = 0;
  {
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      bndry_func(geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchSingleLevel(
      *a_force[lev], amrex::IntVect(nGrowForce), a_time, {a_force[lev]},
      {a_time}, 0, 0, nComp, geom[lev], bndry_func, 0);
  }
  for (lev = 1; lev <= finest_level; ++lev) {
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      crse_bndry_func(
        geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      fine_bndry_func(
        geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    amrex::Interpolater* mapper = &amrex::pc_interp;
    FillPatchTwoLevels(
      *a_force[lev], amrex::IntVect(nGrowForce), a_time, {a_force[lev - 1]},
      {a_time}, {a_force[lev]}, {a_time}, 0, 0, nComp, geom[lev - 1], geom[lev],
      crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1), mapper,
      {m_bcrec_force}, 0);
  }
}

// Fill the gradp
void
PeleLM::fillpatch_gradp(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_gp,
  const int nGhost)
{
  if (lev == 0) {
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      bndry_func(geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchSingleLevel(
      a_gp, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev]->gp), &(m_leveldata_new[lev]->gp)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, AMREX_SPACEDIM, geom[lev], bndry_func,
      0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      crse_bndry_func(
        geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      fine_bndry_func(
        geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchTwoLevels(
      a_gp, amrex::IntVect(nGhost), a_time,
      {&(m_leveldata_old[lev - 1]->gp), &(m_leveldata_new[lev - 1]->gp)},
      {m_t_old[lev - 1], m_t_new[lev - 1]},
      {&(m_leveldata_old[lev]->gp), &(m_leveldata_new[lev]->gp)},
      {m_t_old[lev], m_t_new[lev]}, 0, 0, AMREX_SPACEDIM, geom[lev - 1],
      geom[lev], crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1),
      mapper, {m_bcrec_force}, 0);
  }
}

// Fill the reaction data
void
PeleLM::fillpatch_reaction(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_I_R,
  const int nGhost)
{
  if (lev == 0) {
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      bndry_func(geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchSingleLevel(
      a_I_R, amrex::IntVect(nGhost), a_time, {&(m_leveldatareact[lev]->I_R)},
      {a_time}, 0, 0, nCompIR(), geom[lev], bndry_func, 0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      crse_bndry_func(
        geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      fine_bndry_func(
        geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchTwoLevels(
      a_I_R, amrex::IntVect(nGhost), a_time,
      {&(m_leveldatareact[lev - 1]->I_R)}, {a_time},
      {&(m_leveldatareact[lev]->I_R)}, {a_time}, 0, 0, nCompIR(), geom[lev - 1],
      geom[lev], crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1),
      mapper, {m_bcrec_force}, 0);
  }
}

// Fill functC
void
PeleLM::fillpatch_chemFunctCall(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_fctC,
  const int nGhost)
{
  if (lev == 0) {
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      bndry_func(geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchSingleLevel(
      a_fctC, amrex::IntVect(nGhost), a_time,
      {&(m_leveldatareact[lev]->functC)}, {a_time}, 0, 0, 1, geom[lev],
      bndry_func, 0);
  } else {

    // Interpolator
    auto* mapper = getInterpolator();

    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      crse_bndry_func(
        geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
      fine_bndry_func(
        geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
    FillPatchTwoLevels(
      a_fctC, amrex::IntVect(nGhost), a_time,
      {&(m_leveldatareact[lev - 1]->functC)}, {a_time},
      {&(m_leveldatareact[lev]->functC)}, {a_time}, 0, 0, 1, geom[lev - 1],
      geom[lev], crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1),
      mapper, {m_bcrec_force}, 0);
  }
}

// Fill the state
void
PeleLM::fillcoarsepatch_state(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_state,
  const int nGhost)
{
  AMREX_ASSERT(lev > 0);
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();

  const int nCompState = (m_incompressible) != 0 ? AMREX_SPACEDIM : NVAR;

  fillTurbInflow(a_state, VELX, lev, a_time);

  if (m_use_inlet_from_plane != 0) {
    fillFromRecyclingPlane(a_state, 0, lev);
  }

  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
    crse_bndry_func(
      geom[lev - 1], fetchBCRecArray(0, nCompState),
      PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux,
        static_cast<int>(turb_inflow.is_initialized()),
        static_cast<int>(m_use_inlet_from_plane)});
  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
    fine_bndry_func(
      geom[lev], fetchBCRecArray(0, nCompState),
      PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux,
        static_cast<int>(turb_inflow.is_initialized()),
        static_cast<int>(m_use_inlet_from_plane)});
  InterpFromCoarseLevel(
    a_state, amrex::IntVect(nGhost), a_time, m_leveldata_new[lev - 1]->state, 0,
    0, nCompState, geom[lev - 1], geom[lev], crse_bndry_func, 0,
    fine_bndry_func, 0, refRatio(lev - 1), mapper,
    fetchBCRecArray(0, nCompState), 0);
}

// Fill the auxiliaries
void
PeleLM::fillcoarsepatch_aux(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_aux,
  const int nGhost)
{
  AMREX_ASSERT(lev > 0);
  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();

  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirAux<ProblemSpecificFunctions>>>
    crse_bndry_func(
      geom[lev - 1], fetchBCRecAuxArray(0, m_nAux),
      PeleLMCCFillExtDirAux<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux});
  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirAux<ProblemSpecificFunctions>>>
    fine_bndry_func(
      geom[lev], fetchBCRecAuxArray(0, m_nAux),
      PeleLMCCFillExtDirAux<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux});
  InterpFromCoarseLevel(
    a_aux, amrex::IntVect(nGhost), a_time,
    m_leveldata_new[lev - 1]->auxiliaries, 0, 0, m_nAux, geom[lev - 1],
    geom[lev], crse_bndry_func, 0, fine_bndry_func, 0, refRatio(lev - 1),
    mapper, fetchBCRecAuxArray(0, m_nAux), 0);
}

// Fill the grad P
void
PeleLM::fillcoarsepatch_gradp(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_gp,
  const int nGhost)
{
  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    crse_bndry_func(
      geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    fine_bndry_func(
      geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
  InterpFromCoarseLevel(
    a_gp, amrex::IntVect(nGhost), a_time, m_leveldata_new[lev - 1]->gp, 0, 0,
    AMREX_SPACEDIM, geom[lev - 1], geom[lev], crse_bndry_func, 0,
    fine_bndry_func, 0, refRatio(lev - 1), mapper, {m_bcrec_force}, 0);
}

// Fill the divu
void
PeleLM::fillcoarsepatch_divu(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_divu,
  const int nGhost)
{
  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    crse_bndry_func(
      geom[lev - 1], {m_bcrec_divu}, PeleLMCCFillExtDirDummy{m_nAux});
  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    fine_bndry_func(geom[lev], {m_bcrec_divu}, PeleLMCCFillExtDirDummy{m_nAux});
  InterpFromCoarseLevel(
    a_divu, amrex::IntVect(nGhost), a_time, m_leveldata_new[lev - 1]->divu, 0,
    0, 1, geom[lev - 1], geom[lev], crse_bndry_func, 0, fine_bndry_func, 0,
    refRatio(lev - 1), mapper, {m_bcrec_divu}, 0);
}

// Fill coarse patch of reaction
void
PeleLM::fillcoarsepatch_reaction(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_I_R,
  const int nGhost)
{
  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    crse_bndry_func(
      geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    fine_bndry_func(
      geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
  InterpFromCoarseLevel(
    a_I_R, amrex::IntVect(nGhost), a_time, m_leveldatareact[lev - 1]->I_R, 0, 0,
    nCompIR(), geom[lev - 1], geom[lev], crse_bndry_func, 0, fine_bndry_func, 0,
    refRatio(lev - 1), mapper, {m_bcrec_force}, 0);
}

// Fill coarse patch of chem function call
void
PeleLM::fillcoarsepatch_chemFunctCall(
  const int lev,
  const amrex::Real a_time,
  amrex::MultiFab& a_fctC,
  const int nGhost)
{
  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    crse_bndry_func(
      geom[lev - 1], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<PeleLMCCFillExtDirDummy>>
    fine_bndry_func(
      geom[lev], {m_bcrec_force}, PeleLMCCFillExtDirDummy{m_nAux});
  InterpFromCoarseLevel(
    a_fctC, amrex::IntVect(nGhost), a_time, m_leveldatareact[lev - 1]->functC,
    0, 0, 1, geom[lev - 1], geom[lev], crse_bndry_func, 0, fine_bndry_func, 0,
    refRatio(lev - 1), mapper, {m_bcrec_force}, 0);
}

// Fill the inflow boundary of a velocity MF
// used for velocity projection
void
PeleLM::setInflowBoundaryVel(
  amrex::MultiFab& a_vel, const int lev, const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::setInflowBoundaryVel()");

  const amrex::Real time = getTime(lev, a_time);

  // Create a dummy BCRec from Velocity BCRec keeping only Inflow and set the
  // other to bogus
  auto realVelBCRec = fetchBCRecArray(VELX, AMREX_SPACEDIM);
  amrex::Vector<amrex::BCRec> dummyVelBCRec(AMREX_SPACEDIM);
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    for (int idim2 = 0; idim2 < AMREX_SPACEDIM; ++idim2) {
      if (realVelBCRec[idim].lo(idim2) == amrex::BCType::ext_dir) {
        dummyVelBCRec[idim].setLo(idim2, amrex::BCType::ext_dir);
      } else {
        dummyVelBCRec[idim].setLo(idim2, amrex::BCType::bogus);
      }
      if (realVelBCRec[idim].hi(idim2) == amrex::BCType::ext_dir) {
        dummyVelBCRec[idim].setHi(idim2, amrex::BCType::ext_dir);
      } else {
        dummyVelBCRec[idim].setHi(idim2, amrex::BCType::bogus);
      }
    }
  }

  fillTurbInflow(a_vel, 0, lev, time);

  if (m_use_inlet_from_plane != 0) {
    fillFromRecyclingPlane(a_vel, 0, lev);
  }

  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
    bndry_func(
      geom[lev], dummyVelBCRec,
      PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux,
        static_cast<int>(turb_inflow.is_initialized()),
        static_cast<int>(m_use_inlet_from_plane)});

  bndry_func(a_vel, 0, AMREX_SPACEDIM, a_vel.nGrowVect(), time, 0);

  a_vel.EnforcePeriodicity(geom[lev].periodicity());
}

void
PeleLM::fillTurbInflow(
  amrex::MultiFab& a_vel,
  const int vel_comp,
  const int lev,
  const amrex::Real a_time)
{
  if (turb_inflow.is_initialized()) {

    ProbParm* probparmDD = PeleLM::prob_parm_d;
    ProbParm* probparmDH = PeleLM::prob_parm;

    // Velocity BCs
    auto velBCRec = fetchBCRecArray(VELX, AMREX_SPACEDIM);

    // Copy problem parameter structs to host
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, probparmDD, probparmDD + 1, probparmDH);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(a_vel, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      amrex::Box const& bx = mfi.growntilebox();
      amrex::FArrayBox& data = a_vel[mfi];

      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {

        auto bndryBoxLO =
          amrex::Box(amrex::adjCellLo(geom[lev].Domain(), dir, 4) & bx);
        if (
          velBCRec[0].lo()[dir] == amrex::BCType::ext_dir && bndryBoxLO.ok()) {
          // Create box with ghost cells and set them to zero
          amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
          constexpr int Grow = 4; // Being conservative
          for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            growVect[n] = Grow;
          }
          growVect[dir] = 0;
          amrex::Box modDom = geom[lev].Domain();
          modDom.grow(growVect);
          auto bndryBoxLO_ghost =
            amrex::Box(amrex::adjCellLo(modDom, dir, Grow) & bx);
          data.setVal<amrex::RunOn::Device>(
            0.0, bndryBoxLO_ghost, vel_comp, AMREX_SPACEDIM);

          turb_inflow.add_turb(
            bndryBoxLO, data, 0, geom[lev], a_time, dir,
            amrex::Orientation::low);
        }

        auto bndryBoxHI =
          amrex::Box(amrex::adjCellHi(geom[lev].Domain(), dir, 4) & bx);
        if (
          velBCRec[0].hi()[dir] == amrex::BCType::ext_dir && bndryBoxHI.ok()) {
          // Create box with ghost cells and set them to zero
          amrex::IntVect growVect(amrex::IntVect::TheUnitVector());
          constexpr int Grow = 4;
          for (int n = 0; n < AMREX_SPACEDIM; ++n) {
            growVect[n] = Grow;
          }
          growVect[dir] = 0;
          amrex::Box modDom = geom[lev].Domain();
          modDom.grow(growVect);
          auto bndryBoxHI_ghost =
            amrex::Box(amrex::adjCellHi(modDom, dir, Grow) & bx);
          data.setVal<amrex::RunOn::Device>(
            0.0, bndryBoxHI_ghost, vel_comp, AMREX_SPACEDIM);

          turb_inflow.add_turb(
            bndryBoxHI, data, 0, geom[lev], a_time, dir,
            amrex::Orientation::high);
        }
      }
    }

    // Copy problem parameter structs back to device
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, probparmDH, probparmDH + 1, probparmDD);
  }
}

int
PeleLM::computeRecyclingSrcIndex(int lev) const
{
  const int dir = m_inlet_plane_dir;
  return static_cast<int>(std::lround(
    (m_inlet_plane_position - geom[lev].ProbLo()[dir]) /
      geom[lev].CellSize()[dir] -
    0.5));
}

void
PeleLM::buildRecyclingPlaneStorage()
{
  if (m_use_inlet_from_plane == 0) {
    return;
  }

  const int planeDir = m_inlet_plane_dir;
  const int nlevels = finest_level + 1;

  // The slab BoxArray at each level depends only on the (fixed) domain,
  // maxGridSize, and srcIndex, so it is invariant across regrids. The level's
  // grids and the DistributionMapping may change, but the running mean is a
  // physical-space quantity and can be carried through with a ParallelCopy.
  // Stash the existing means before reallocating so we can re-deposit them
  // into the new MultiFabs.
  auto saved_mean = std::move(m_inlet_recycling.mean_src);
  const bool saved_initialized = m_inlet_recycling.initialized;
  const int saved_n_samples = m_inlet_recycling.n_samples;

  m_inlet_recycling.u_src.clear();
  m_inlet_recycling.fluct_src.clear();
  m_inlet_recycling.mean_src.clear();
  m_inlet_recycling.u_src.resize(nlevels);
  m_inlet_recycling.mean_src.resize(nlevels);
  m_inlet_recycling.fluct_src.resize(nlevels);
#ifdef AMREX_USE_EB
  m_inlet_recycling.mask.clear();
  m_inlet_recycling.mask.resize(nlevels);
#endif

  // If a level was added by this regrid, we have no saved mean for it. Falling
  // back to a fresh reseed is safer than zero-mean (which would inject the
  // full instantaneous velocity, not a fluctuation, on the new level).
  bool full_reseed = false;

  for (int lev = 0; lev < nlevels; ++lev) {
    const int srcIndex = computeRecyclingSrcIndex(lev);
    const amrex::Box& domain = geom[lev].Domain();
    AMREX_ALWAYS_ASSERT(
      srcIndex >= domain.smallEnd(planeDir) &&
      srcIndex <= domain.bigEnd(planeDir));

    // Thin slab spanning the entire transverse cross-section at srcIndex.
    // Build it by intersecting the slab with this level's existing grids so
    // the slab MultiFabs inherit the same processor ownership as the source
    // state data and subsequent ParallelCopy operations preserve locality.
    amrex::Box slab = domain;
    slab.setSmall(planeDir, srcIndex);
    slab.setBig(planeDir, srcIndex);

    const amrex::BoxArray& level_ba = boxArray(lev);
    const amrex::DistributionMapping& level_dm = DistributionMap(lev);
    const auto& level_pmap = level_dm.ProcessorMap();
    amrex::BoxList slab_bl;
    amrex::Vector<int> slab_pmap;
    for (int ibox = 0; ibox < level_ba.size(); ++ibox) {
      amrex::Box isect = level_ba[ibox] & slab;
      if (isect.ok()) {
        slab_bl.push_back(isect);
        slab_pmap.push_back(level_pmap[ibox]);
      }
    }
    AMREX_ALWAYS_ASSERT(!slab_pmap.empty());
    amrex::BoxArray slab_ba(slab_bl);
    amrex::DistributionMapping slab_dm(slab_pmap);

    m_inlet_recycling.u_src[lev] =
      std::make_unique<amrex::MultiFab>(slab_ba, slab_dm, AMREX_SPACEDIM, 0);
    m_inlet_recycling.mean_src[lev] =
      std::make_unique<amrex::MultiFab>(slab_ba, slab_dm, AMREX_SPACEDIM, 0);
    m_inlet_recycling.fluct_src[lev] =
      std::make_unique<amrex::MultiFab>(slab_ba, slab_dm, AMREX_SPACEDIM, 0);

    // No fluctuation until the next snapshot has populated u_src and
    // recomputed it.
    m_inlet_recycling.fluct_src[lev]->setVal(0.0);

#ifdef AMREX_USE_EB
    // Build a fresh EB factory at the slab's BoxArray to obtain per-cell
    // flags at this level's resolution, then translate to a 0/1 mask:
    // 0 = EB-covered (excluded from the running mean; zero fluctuation),
    // 1 = regular or cut (included).
    auto slab_eb_factory = amrex::makeEBFabFactory(
      geom[lev], slab_ba, slab_dm, {AMREX_D_DECL(0, 0, 0)},
      amrex::EBSupport::basic);
    const auto& flags = slab_eb_factory->getMultiEBCellFlagFab();

    m_inlet_recycling.mask[lev] =
      std::make_unique<amrex::iMultiFab>(slab_ba, slab_dm, 1, 0);
    auto& mask_lev = *m_inlet_recycling.mask[lev];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mask_lev, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      const auto& flagarr = flags.const_array(mfi);
      auto const& mask_arr = mask_lev.array(mfi);
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          mask_arr(i, j, k) = flagarr(i, j, k).isCovered() ? 0 : 1;
        });
    }
#endif

    if (
      saved_initialized && lev < static_cast<int>(saved_mean.size()) &&
      saved_mean[lev]) {
      m_inlet_recycling.mean_src[lev]->ParallelCopy(
        *saved_mean[lev], 0, 0, AMREX_SPACEDIM);
    } else {
      // Either we never had a mean, or this level didn't exist before.
      m_inlet_recycling.mean_src[lev]->setVal(0.0);
      if (saved_initialized) {
        full_reseed = true;
      }
    }
  }

  if (saved_initialized && !full_reseed) {
    // Carry the running statistics forward.
    m_inlet_recycling.initialized = true;
    m_inlet_recycling.n_samples = saved_n_samples;
  } else {
    // No prior mean, or a new level appeared and we don't have a mean for it
    // — let the next snapshot reseed cleanly.
    m_inlet_recycling.initialized = false;
    m_inlet_recycling.n_samples = 0;
  }
}

void
PeleLM::updateRecyclingPlaneSnapshot()
{
  if (m_use_inlet_from_plane == 0) {
    return;
  }

  if (m_recyclingNeedsRebuild != 0) {
    buildRecyclingPlaneStorage();
    m_recyclingNeedsRebuild = 0;
  }

  // The new-time velocity is what we sample; ensure storage exists.
  AMREX_ASSERT(
    static_cast<int>(m_inlet_recycling.u_src.size()) == finest_level + 1);

  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  auto velBCRec = fetchBCRecArray(VELX, AMREX_SPACEDIM);
  // The time used by InterpFromCoarseLevel is mostly informational here
  // (the slab is in the interior, so PhysBCFunct calls on its temporaries
  // are no-ops in practice); pass the new time for consistency.
  const amrex::Real a_time = m_cur_time;

  for (int lev = 0; lev <= finest_level; ++lev) {
    auto& u_src = *m_inlet_recycling.u_src[lev];
    const auto& state_lev = m_leveldata_new[lev]->state;

    // Every level's state covers some, but in general not all, of the
    // transverse cross-section at srcIndex. Strategy:
    //   1. Lev > 0: interpolate from the (already-filled) coarser slab to
    //      get a complete coverage at this level's resolution.
    //   2. Overwrite from this level's own state where it covers, using the
    //      higher-resolution data.
    //   NOTE: Relies on data filled at next coarse level, filled by this loop.
    if (lev > 0) {
      // Use piecewise-constant interpolation: the source slab is one cell
      // thick along planeDir, so any stencil-based interpolator (e.g.,
      // cell-conservative linear) would read garbage from the slab's
      // planeDir ghost cells. PCInterp has no transverse stencil and is
      // adequate for injecting a fluctuation field across a refinement
      // boundary.
      // NOTE: Declared as InterpBase* (not auto* / MFPCInterp*) so AMReX's
      // FillPatchInterp dispatches through its runtime dynamic_cast path
      // and picks the MultiFab-based interpolater entry point.
      amrex::InterpBase* mapper = &amrex::mf_pc_interp;
      amrex::PhysBCFunct<amrex::GpuBndryFuncFab<
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
        crse_bndry_func(
          geom[lev - 1], velBCRec,
          PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
            lprobparm, lpmfdata, m_nAux,
            static_cast<int>(turb_inflow.is_initialized()),
            static_cast<int>(m_use_inlet_from_plane)});
      amrex::PhysBCFunct<amrex::GpuBndryFuncFab<
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
        fine_bndry_func(
          geom[lev], velBCRec,
          PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
            lprobparm, lpmfdata, m_nAux,
            static_cast<int>(turb_inflow.is_initialized()),
            static_cast<int>(m_use_inlet_from_plane)});
      amrex::InterpFromCoarseLevel(
        u_src, amrex::IntVect(0), a_time, *m_inlet_recycling.u_src[lev - 1], 0,
        0, AMREX_SPACEDIM, geom[lev - 1], geom[lev], crse_bndry_func, 0,
        fine_bndry_func, 0, refRatio(lev - 1), mapper, velBCRec, 0);
    }
    // Same-level data overrides the coarse-interpolated baseline anywhere
    // this level's grids cover the slab.
    u_src.ParallelCopy(
      state_lev, VELX, 0, AMREX_SPACEDIM, 0, 0, geom[lev].periodicity());

#ifdef AMREX_USE_EB
    // EB-covered cells contain undefined storage; zero them in u_src so they
    // contribute nothing to the running mean (since 0 is the steady value
    // there) and produce a zero fluctuation downstream.
    if (
      lev < static_cast<int>(m_inlet_recycling.mask.size()) &&
      m_inlet_recycling.mask[lev]) {
      const auto& mask_lev = *m_inlet_recycling.mask[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(u_src, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        const auto& mask_arr = mask_lev.const_array(mfi);
        auto const& u_arr = u_src.array(mfi);
        amrex::ParallelFor(
          bx, AMREX_SPACEDIM,
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (mask_arr(i, j, k) == 0) {
              u_arr(i, j, k, n) = 0.0;
            }
          });
      }
    }
#endif
  }

  // Update the running mean and store the current fluctuation.
  ++m_inlet_recycling.n_samples;

  if (!m_inlet_recycling.initialized) {
    // Seed: <u> = u_0; fluctuation defined as zero on the seeding sample.
    for (int lev = 0; lev <= finest_level; ++lev) {
      amrex::MultiFab::Copy(
        *m_inlet_recycling.mean_src[lev], *m_inlet_recycling.u_src[lev], 0, 0,
        AMREX_SPACEDIM, 0);
      m_inlet_recycling.fluct_src[lev]->setVal(0.0);
    }
    m_inlet_recycling.initialized = true;
    return;
  }

  // alpha for exponential moving average; if no window is set, fall back to a
  // cumulative average via 1/n.
  amrex::Real alpha;
  if (m_inlet_plane_avg_window > 0.0) {
    alpha = std::min(amrex::Real(1.0), m_dt / m_inlet_plane_avg_window);
    if (alpha == amrex::Real(1.0)) {
      amrex::Print() << "WARNING: Clipped recycle averaging window will give "
                        "no fluctuations.\n";
    }
  } else {
    alpha =
      1.0 / static_cast<amrex::Real>(std::max(1, m_inlet_recycling.n_samples));
  }
  const amrex::Real one_minus_alpha = 1.0 - alpha;

  for (int lev = 0; lev <= finest_level; ++lev) {
    auto& mean = *m_inlet_recycling.mean_src[lev];
    auto& u_src = *m_inlet_recycling.u_src[lev];
    auto& fluct = *m_inlet_recycling.fluct_src[lev];

    // mean = (1 - alpha) * mean + alpha * u_src
    amrex::MultiFab::LinComb(
      mean, one_minus_alpha, mean, 0, alpha, u_src, 0, 0, AMREX_SPACEDIM, 0);
    // fluct = u_src - mean
    amrex::MultiFab::LinComb(
      fluct, 1.0, u_src, 0, -1.0, mean, 0, 0, AMREX_SPACEDIM, 0);
  }
}

void
PeleLM::fillFromRecyclingPlane(amrex::MultiFab& a_vel, int vel_comp, int lev)
{
  // NOTE: Fluctuation data injected is refreshed only once per time step
  // (rather than per SDC iteration)

  if (m_use_inlet_from_plane == 0) {
    return;
  }

  // Recycling currently only handles cell-centered velocity data; staggered
  // (e.g., MAC) layouts would need a different shift convention.
  AMREX_ASSERT(a_vel.boxArray().ixType().cellCentered());

  // Storage may not yet exist (first call before updateRecyclingPlaneSnapshot,
  // or this level didn't exist when the snapshot last ran). Fall back to the
  // standard ext_dir fill silently.
  if (
    lev >= static_cast<int>(m_inlet_recycling.fluct_src.size()) ||
    !m_inlet_recycling.fluct_src[lev]) {
    return;
  }
  // Don't inject anything until the running mean has had a chance to settle.
  if (
    !m_inlet_recycling.initialized ||
    m_inlet_recycling.n_samples <= m_inlet_plane_warmup_steps) {
    return;
  }

  const int planeDir = m_inlet_plane_dir;
  const int srcIndex = computeRecyclingSrcIndex(lev);
  const amrex::Box& domain = geom[lev].Domain();
  if (
    srcIndex < domain.smallEnd(planeDir) ||
    srcIndex > domain.bigEnd(planeDir)) {
    amrex::Print() << "[fillFromRecyclingPlane] lev " << lev
                   << " plane outside domain\n";
    return;
  }

  auto velBCRec = fetchBCRecArray(VELX, AMREX_SPACEDIM);
  const amrex::BoxArray& ba = a_vel.boxArray();
  const int nGrowDest = a_vel.nGrow();

  // Recycling injects all AMREX_SPACEDIM velocity components as a single
  // operation, so it only makes sense when every velocity component on the
  // affected face is ext_dir. Disagreement is almost certainly a malformed
  // input, so abort rather than silently picking a side.
  auto faceIsExtDir = [&](amrex::Orientation::Side side) {
    int n_extdir = 0;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      const int bc = (side == amrex::Orientation::low)
                       ? velBCRec[idim].lo()[planeDir]
                       : velBCRec[idim].hi()[planeDir];
      if (bc == amrex::BCType::ext_dir) {
        ++n_extdir;
      }
    }
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      n_extdir == 0 || n_extdir == AMREX_SPACEDIM,
      "fillFromRecyclingPlane: velocity components disagree on ext_dir for "
      "the recycling face; all components must share the same BC.");
    return n_extdir == AMREX_SPACEDIM;
  };

  // Decide which sides of planeDir need recycling. Uses only global
  // quantities (BCRec + domain + ba), so the result is identical on
  // every rank and the subsequent ParallelAdd calls remain collective.
  bool need_lo = false;
  bool need_hi = false;
  amrex::BoxArray grown_ba(ba);
  grown_ba.grow(nGrowDest);
  if (faceIsExtDir(amrex::Orientation::low)) {
    need_lo =
      grown_ba.intersects(amrex::adjCellLo(domain, planeDir, nGrowDest));
  }
  if (faceIsExtDir(amrex::Orientation::high)) {
    need_hi =
      grown_ba.intersects(amrex::adjCellHi(domain, planeDir, nGrowDest));
  }

  if (!need_lo && !need_hi) {
    return;
  }

  // Shift the cached fluctuation MultiFab onto each inflow ghost layer and
  // COPY it to the destination. The standard ext_dir fill will add
  // the inlet mean profile; this only contributes the zero-mean fluctuation.
  amrex::MultiFab& fluct = *m_inlet_recycling.fluct_src[lev];

  auto copy_shifted_fluct = [&](const amrex::IntVect& shift) {
    amrex::BoxArray shifted_ba(fluct.boxArray());
    shifted_ba.shift(shift);

    amrex::MultiFab shifted_fluct(
      shifted_ba, fluct.DistributionMap(), fluct.nComp(), 0, amrex::MFInfo(),
      fluct.Factory());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(fluct); mfi.isValid(); ++mfi) {
      const amrex::Box& src_bx = fluct[mfi].box();
      const amrex::Box dst_bx = amrex::shift(src_bx, shift);
      shifted_fluct[mfi].copy(fluct[mfi], src_bx, 0, dst_bx, 0, fluct.nComp());
    }

    a_vel.ParallelCopy(
      shifted_fluct, 0, vel_comp, AMREX_SPACEDIM, 0, nGrowDest);
  };

  if (need_lo) {
    for (int g = 1; g <= nGrowDest; ++g) {
      const int nshift = srcIndex - domain.smallEnd(planeDir) + g;
      const auto shift = amrex::BASISV(planeDir) * nshift;
      copy_shifted_fluct(shift);
    }
  }
  if (need_hi) {
    for (int g = 1; g <= nGrowDest; ++g) {
      const int nshift = domain.bigEnd(planeDir) - srcIndex + g;
      const auto shift = amrex::BASISV(planeDir) * nshift;
      copy_shifted_fluct(shift);
    }
  }
}
