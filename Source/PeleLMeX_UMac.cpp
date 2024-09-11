#include <PeleLMeX.H>
#include <PeleLMeX_Utils.H>
#include <hydro_utils.H>
#include <AMReX_FillPatchUtil.H>
#include <PeleLMeX_BCfill.H>

using namespace amrex;

void
PeleLM::predictVelocity(std::unique_ptr<AdvanceAdvData>& advData)
{
  BL_PROFILE("PeleLMeX::predictVelocity()");

  // set umac boundaries to zero
  if (advData->umac[0][0].nGrow() > 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      AMREX_D_TERM(advData->umac[lev][0].setBndry(0.0);
                   , advData->umac[lev][1].setBndry(0.0);
                   , advData->umac[lev][2].setBndry(0.0););
    }
  }

  //----------------------------------------------------------------
  // Get viscous forces
  int nGrow_force = 1;
  Vector<MultiFab> divtau(finest_level + 1);
  Vector<MultiFab> velForces(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    divtau[lev].define(
      grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), Factory(lev));
    velForces[lev].define(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force, MFInfo(),
      Factory(lev));
  }
  int use_density = 0;
  computeDivTau(AmrOldTime, GetVecOfPtrs(divtau), use_density);

  //----------------------------------------------------------------
  // Gather all the velocity forces
  // F = [ (gravity+...) - gradP + divTau ] / rho
  int add_gradP = 1;
  getVelForces(
    AmrOldTime, GetVecOfPtrs(divtau), GetVecOfPtrs(velForces), nGrow_force,
    add_gradP);

  //----------------------------------------------------------------
  // Predict face velocities at t^{n+1/2} with Godunov
  auto bcRecVel = fetchBCRecArray(VELX, AMREX_SPACEDIM);
  auto bcRecVel_d = convertToDeviceVector(bcRecVel);
  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, AmrOldTime);

    //----------------------------------------------------------------
#ifdef AMREX_USE_EB
    const auto& ebfact = EBFactory(lev);
#endif

    HydroUtils::ExtrapVelToFaces(
      ldata_p->state, velForces[lev],
      AMREX_D_DECL(
        advData->umac[lev][0], advData->umac[lev][1], advData->umac[lev][2]),
      bcRecVel, bcRecVel_d.dataPtr(), geom[lev], m_dt,
#ifdef AMREX_USE_EB
      ebfact,
#endif
      m_Godunov_ppm != 0, m_Godunov_ForceInTrans != 0, m_predict_advection_type,
      m_Godunov_ppm_limiter);
  }
}

void
PeleLM::createMACRHS(std::unique_ptr<AdvanceAdvData>& advData)
{
  BL_PROFILE("PeleLMeX::createMACRHS()");

  for (int lev = 0; lev <= finest_level; ++lev) {
    Real halftime = 0.5 * (m_t_old[lev] + m_t_new[lev]);
    fillpatch_divu(lev, halftime, advData->mac_divu[lev], m_nGrowAdv);
  }
}

void
PeleLM::addChiIncrement(
  int a_sdcIter,
  const TimeStamp& a_time,
  std::unique_ptr<AdvanceAdvData>& advData)
{
  BL_PROFILE("PeleLMeX::addChiIncrement()");

  int nGrow = m_nGrowAdv;
  Vector<MultiFab> chiIncr(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    chiIncr[lev].define(
      grids[lev], dmap[lev], 1, nGrow, MFInfo(), Factory(lev));
  }

  // Update the thermodynamic pressure
  setThermoPress(a_time);

  // Compute the pressure drift term
  calc_dPdt(a_time, GetVecOfPtrs(chiIncr));

  // Add chiIncr to chi and add chi to mac_divu
  // Both mac_divu and chiIncr have properly filled ghost cells -> work on
  // grownbox
  for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(advData->chi[lev], TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const Box& gbx = mfi.growntilebox();
      auto const& chiInc_ar = chiIncr[lev].const_array(mfi);
      auto const& chi_ar = advData->chi[lev].array(mfi);
      auto const& mac_divu_ar = advData->mac_divu[lev].array(mfi);
      if (m_chi_correction_type == ChiCorrectionType::DivuFirstIter) {
        amrex::ParallelFor(
          gbx, [chi_ar, chiInc_ar, mac_divu_ar,
                a_sdcIter] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (a_sdcIter == 1) {
              chi_ar(i, j, k) = chiInc_ar(i, j, k) + mac_divu_ar(i, j, k);
            } else {
              chi_ar(i, j, k) += chiInc_ar(i, j, k);
            }
            mac_divu_ar(i, j, k) = chi_ar(i, j, k);
          });
      } else if (m_chi_correction_type == ChiCorrectionType::NoDivu) {
        amrex::ParallelFor(
          gbx, [chi_ar, chiInc_ar, mac_divu_ar,
                a_sdcIter] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (a_sdcIter == 1) {
              chi_ar(i, j, k) = chiInc_ar(i, j, k);
            } else {
              chi_ar(i, j, k) += chiInc_ar(i, j, k);
            }
            mac_divu_ar(i, j, k) = chi_ar(i, j, k);
          });
      } else { // Default: use updated divu every iteration
        amrex::ParallelFor(
          gbx, [chi_ar, chiInc_ar, mac_divu_ar,
                a_sdcIter] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (a_sdcIter == 1) {
              chi_ar(i, j, k) = chiInc_ar(i, j, k);
            } else {
              chi_ar(i, j, k) += chiInc_ar(i, j, k);
            }
            mac_divu_ar(i, j, k) += chi_ar(i, j, k);
          });
      }
    }
  }
  if (m_print_chi_convergence) {
    amrex::Real max_corr =
      MLNorm0(GetVecOfConstPtrs(chiIncr)) * m_dt / m_dpdtFactor;
    amrex::Print() << "      Before SDC " << a_sdcIter
                   << ": max relative P mismatch is " << max_corr << std::endl;
  }
}

void
PeleLM::macProject(
  const TimeStamp& a_time,
  std::unique_ptr<AdvanceAdvData>& advData,
  const Vector<MultiFab*>& a_divu)
{
  BL_PROFILE("PeleLMeX::macProject()");

  int has_divu = static_cast<int>(!a_divu.empty());

  // Get face rho inv
  auto bcRec = fetchBCRecArray(DENSITY, 1);
  Vector<Array<MultiFab, AMREX_SPACEDIM>> rho_inv(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (m_incompressible != 0) {
      Real rhoInv = m_dt / (2.0 * m_rho);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        rho_inv[lev][idim].define(
          amrex::convert(grids[lev], IntVect::TheDimensionVector(idim)),
          dmap[lev], 1, 0, MFInfo(), Factory(lev));
        rho_inv[lev][idim].setVal(rhoInv);
      }
    } else {
      auto* ldata_p = getLevelDataPtr(lev, a_time);
      int doZeroVisc = 0;
      rho_inv[lev] =
        getDiffusivity(lev, DENSITY, 1, doZeroVisc, {bcRec}, ldata_p->state);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        rho_inv[lev][idim].invert(m_dt / 2.0, 0);
        rho_inv[lev][idim].FillBoundary(geom[lev].periodicity());
      }
    }
  }

  // For closed chamber, compute change in chamber pressure
  Real Sbar = 0.0;
  if ((m_closed_chamber != 0) && (m_incompressible == 0)) {
    Sbar = adjustPandDivU(advData);
  }

  if (macproj->needInitialization()) {
    LPInfo lpInfo;
    lpInfo.setMaxCoarseningLevel(m_mac_mg_max_coarsening_level);
    macproj->initProjector(lpInfo, GetVecOfArrOfConstPtrs(rho_inv));
    macproj->setDomainBC(
      getMACProjectionBC(Orientation::low),
      getMACProjectionBC(Orientation::high));
#ifdef AMREX_USE_HYPRE
    macproj->getMLMG().setHypreOptionsNamespace(m_hypre_namespace_mac);
#endif
  } else {
    macproj->updateBeta(GetVecOfArrOfConstPtrs(rho_inv));
  }

  // set MAC velocity and projection RHS
  macproj->getLinOp().setMaxOrder(m_mac_max_order);
  macproj->setUMAC(GetVecOfArrOfPtrs(advData->umac));
  if (has_divu != 0) {
    macproj->setDivU(GetVecOfConstPtrs(a_divu));
  }

  // Project
  macproj->project(m_mac_mg_rtol, m_mac_mg_atol);

  // Restore mac_divu
  if ((m_closed_chamber != 0) && (m_incompressible == 0)) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      a_divu[lev]->plus(Sbar, 0, 1);
    }
  }

  // FillBoundary umac
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (lev > 0) {
      // We need to fill the MAC velocities outside the fine region so we can
      // use them in the Godunov method
      IntVect rr = geom[lev].Domain().size() / geom[lev - 1].Domain().size();
      create_constrained_umac_grown(
        m_nGrowMAC, &geom[lev - 1], &geom[lev],
        GetArrOfPtrs(advData->umac[lev - 1]), GetArrOfPtrs(advData->umac[lev]),
        rr);
    } else {
      AMREX_D_TERM(
        advData->umac[lev][0].FillBoundary(geom[lev].periodicity());
        , advData->umac[lev][1].FillBoundary(geom[lev].periodicity());
        , advData->umac[lev][2].FillBoundary(geom[lev].periodicity()));
    }
  }
}

void
PeleLM::create_constrained_umac_grown(
  int a_nGrow,
  const Geometry* crse_geom,
  const Geometry* fine_geom,
  Array<MultiFab*, AMREX_SPACEDIM> u_mac_crse,
  Array<MultiFab*, AMREX_SPACEDIM> u_mac_fine,
  const IntVect& crse_ratio)
{
  // Divergence preserving interp
  Interpolater* mapper = &face_divfree_interp;

  // Set BCRec for Umac
  Vector<BCRec> bcrec(1);
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    if (crse_geom->isPeriodic(idim)) {
      bcrec[0].setLo(idim, BCType::int_dir);
      bcrec[0].setHi(idim, BCType::int_dir);
    } else {
      bcrec[0].setLo(idim, BCType::foextrap);
      bcrec[0].setHi(idim, BCType::foextrap);
    }
  }
  Array<Vector<BCRec>, AMREX_SPACEDIM> bcrecArr = {
    AMREX_D_DECL(bcrec, bcrec, bcrec)};

  PhysBCFunct<GpuBndryFuncFab<umacFill>> crse_bndry_func(
    *crse_geom, bcrec, umacFill{});
  Array<PhysBCFunct<GpuBndryFuncFab<umacFill>>, AMREX_SPACEDIM> cbndyFuncArr = {
    AMREX_D_DECL(crse_bndry_func, crse_bndry_func, crse_bndry_func)};

  PhysBCFunct<GpuBndryFuncFab<umacFill>> fine_bndry_func(
    *fine_geom, bcrec, umacFill{});
  Array<PhysBCFunct<GpuBndryFuncFab<umacFill>>, AMREX_SPACEDIM> fbndyFuncArr = {
    AMREX_D_DECL(fine_bndry_func, fine_bndry_func, fine_bndry_func)};

  // Use piecewise constant interpolation in time, so create dummy variable for
  // time
  Real dummy = 0.;
  FillPatchTwoLevels(
    u_mac_fine, IntVect(a_nGrow), dummy, {u_mac_crse}, {dummy}, {u_mac_fine},
    {dummy}, 0, 0, 1, *crse_geom, *fine_geom, cbndyFuncArr, 0, fbndyFuncArr, 0,
    crse_ratio, mapper, bcrecArr, 0);
}

Array<LinOpBCType, AMREX_SPACEDIM>
PeleLM::getMACProjectionBC(Orientation::Side a_side)
{

  Array<LinOpBCType, AMREX_SPACEDIM> r;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (Geom(0).isPeriodic(idim)) {
      r[idim] = LinOpBCType::Periodic;
    } else {
      auto physbc =
        (a_side == Orientation::low) ? m_phys_bc.lo(idim) : m_phys_bc.hi(idim);
      if (physbc == amrex::PhysBCType::outflow) {
        r[idim] = LinOpBCType::Dirichlet;
      } else {
        r[idim] = LinOpBCType::Neumann;
      }
    }
  }
  return r;
}
