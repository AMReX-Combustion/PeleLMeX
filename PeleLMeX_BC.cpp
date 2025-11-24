#include <PeleLMeX.H>
#include <PeleLMeX_BCfill.H>
#include <AMReX_FillPatchUtil.H>
#include <memory>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
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
  constexpr int sizeForceBC = amrex::max(AMREX_SPACEDIM, NUM_SPECIES + 2);
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

  if (lev == 0) {
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
      bndry_func(
        geom[lev], fetchBCRecArray(0, nCompState),
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux,
          static_cast<int>(turb_inflow.is_initialized())});
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
          static_cast<int>(turb_inflow.is_initialized())});
    amrex::PhysBCFunct<
      amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
      fine_bndry_func(
        geom[lev], fetchBCRecArray(0, nCompState),
        PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
          lprobparm, lpmfdata, m_nAux,
          static_cast<int>(turb_inflow.is_initialized())});
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

  // Interpolator
  auto* mapper = getInterpolator(m_regrid_interp_method);

  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
    crse_bndry_func(
      geom[lev - 1], fetchBCRecArray(0, nCompState),
      PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux,
        static_cast<int>(turb_inflow.is_initialized())});
  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
    fine_bndry_func(
      geom[lev], fetchBCRecArray(0, nCompState),
      PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux,
        static_cast<int>(turb_inflow.is_initialized())});
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

  ProbParm const* lprobparm = prob_parm_d;
  auto const* lpmfdata = pmf_data.device_parm();
  amrex::PhysBCFunct<
    amrex::GpuBndryFuncFab<PeleLMCCFillExtDirState<ProblemSpecificFunctions>>>
    bndry_func(
      geom[lev], dummyVelBCRec,
      PeleLMCCFillExtDirState<ProblemSpecificFunctions>{
        lprobparm, lpmfdata, m_nAux,
        static_cast<int>(turb_inflow.is_initialized())});

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
