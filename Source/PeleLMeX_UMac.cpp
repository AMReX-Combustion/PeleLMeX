#include <PeleLMeX.H>
#include <PeleLMeX_Utils.H>
#include <hydro_utils.H>
#include <AMReX_FillPatchUtil.H>
#include <PeleLMeX_BCfill.H>

void
PeleLM::predictVelocity(const std::unique_ptr<AdvanceAdvData>& advData)
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
  constexpr int nGrow_force = 1;
  amrex::Vector<amrex::MultiFab> divtau;
  divtau.reserve(finest_level + 1);
  amrex::Vector<amrex::MultiFab> velForces;
  velForces.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    divtau.emplace_back(
      grids[lev], dmap[lev], AMREX_SPACEDIM, 0, amrex::MFInfo(), Factory(lev));
    velForces.emplace_back(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force, amrex::MFInfo(),
      Factory(lev));
  }
  constexpr int use_density = 0;
  computeDivTau(AmrOldTime, GetVecOfPtrs(divtau), use_density);

  //----------------------------------------------------------------
  // Gather all the velocity forces
  // F = [ (gravity+...) - gradP + divTau ] / rho
  constexpr int add_gradP = 1;
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
      (m_useEBinflow != 0)
        ? getEBState(lev, VELX, AMREX_SPACEDIM, AmrOldTime).get()
        : nullptr,
#endif
      m_Godunov_ppm != 0, m_Godunov_ForceInTrans != 0, m_predict_advection_type,
      m_Godunov_ppm_limiter);
  }
}

void
PeleLM::createMACRHS(const std::unique_ptr<AdvanceAdvData>& advData)
{
  BL_PROFILE("PeleLMeX::createMACRHS()");

  for (int lev = 0; lev <= finest_level; ++lev) {
    amrex::Real halftime = 0.5 * (m_t_old[lev] + m_t_new[lev]);
    fillpatch_divu(lev, halftime, advData->mac_divu[lev], m_nGrowAdv);
  }
}

void
PeleLM::addChiIncrement(
  const int a_sdcIter,
  const TimeStamp a_time,
  const std::unique_ptr<AdvanceAdvData>& advData)
{
  BL_PROFILE("PeleLMeX::addChiIncrement()");

  const int nGrow = m_nGrowAdv;
  amrex::Vector<amrex::MultiFab> chiIncr;
  chiIncr.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    chiIncr.emplace_back(
      grids[lev], dmap[lev], 1, nGrow, amrex::MFInfo(), Factory(lev));
  }

  // Update the thermodynamic pressure
  setThermoPress(a_time);

  // Compute the pressure drift term
  calc_dPdt(a_time, GetVecOfPtrs(chiIncr));

  // Add chiIncr to chi and add chi to mac_divu
  // Both mac_divu and chiIncr have properly filled ghost cells -> work on
  // grownbox

  if (a_sdcIter == 1) {
    switch (m_chi_correction_type) {
    case ChiCorrectionType::DivuFirstIter: {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& chiInc_ma = chiIncr[lev].const_arrays();
        auto const& chi_ma = advData->chi[lev].arrays();
        auto const& mac_divu_ma = advData->mac_divu[lev].arrays();
        amrex::ParallelFor(
          advData->chi[lev], advData->chi[lev].nGrowVect(),
          [chi_ma, chiInc_ma, mac_divu_ma] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            chi_ma[box_no](i, j, k) =
              chiInc_ma[box_no](i, j, k) + mac_divu_ma[box_no](i, j, k);
            mac_divu_ma[box_no](i, j, k) = chi_ma[box_no](i, j, k);
          });
        // Shift outside?
        amrex::Gpu::streamSynchronize();
      }
      break;
    }
    case ChiCorrectionType::NoDivu: {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& chiInc_ma = chiIncr[lev].const_arrays();
        auto const& chi_ma = advData->chi[lev].arrays();
        auto const& mac_divu_ma = advData->mac_divu[lev].arrays();
        amrex::ParallelFor(
          advData->chi[lev], advData->chi[lev].nGrowVect(),
          [chi_ma, chiInc_ma, mac_divu_ma] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            chi_ma[box_no](i, j, k) = chiInc_ma[box_no](i, j, k);
            mac_divu_ma[box_no](i, j, k) = chi_ma[box_no](i, j, k);
          });
        // Shift outside?
        amrex::Gpu::streamSynchronize();
      }
      break;
    }
    default: {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& chiInc_ma = chiIncr[lev].const_arrays();
        auto const& chi_ma = advData->chi[lev].arrays();
        auto const& mac_divu_ma = advData->mac_divu[lev].arrays();
        amrex::ParallelFor(
          advData->chi[lev], advData->chi[lev].nGrowVect(),
          [chi_ma, chiInc_ma, mac_divu_ma] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            chi_ma[box_no](i, j, k) = chiInc_ma[box_no](i, j, k);
            mac_divu_ma[box_no](i, j, k) += chi_ma[box_no](i, j, k);
          });
        // Shift outside?
        amrex::Gpu::streamSynchronize();
      }
    }
    }
  } else {
    if (
      m_chi_correction_type == ChiCorrectionType::DivuFirstIter ||
      m_chi_correction_type == ChiCorrectionType::NoDivu) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& chiInc_ma = chiIncr[lev].const_arrays();
        auto const& chi_ma = advData->chi[lev].arrays();
        auto const& mac_divu_ma = advData->mac_divu[lev].arrays();
        amrex::ParallelFor(
          advData->chi[lev], advData->chi[lev].nGrowVect(),
          [chi_ma, chiInc_ma, mac_divu_ma] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            chi_ma[box_no](i, j, k) += chiInc_ma[box_no](i, j, k);
            mac_divu_ma[box_no](i, j, k) = chi_ma[box_no](i, j, k);
          });
        // Shift outside?
        amrex::Gpu::streamSynchronize();
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto const& chiInc_ma = chiIncr[lev].const_arrays();
        auto const& chi_ma = advData->chi[lev].arrays();
        auto const& mac_divu_ma = advData->mac_divu[lev].arrays();
        amrex::ParallelFor(
          advData->chi[lev], advData->chi[lev].nGrowVect(),
          [chi_ma, chiInc_ma, mac_divu_ma] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            chi_ma[box_no](i, j, k) += chiInc_ma[box_no](i, j, k);
            mac_divu_ma[box_no](i, j, k) += chi_ma[box_no](i, j, k);
          });
        // Shift outside?
        amrex::Gpu::streamSynchronize();
      }
    }
  }

  if (m_print_chi_convergence) {
    const amrex::Real max_corr =
      MLNorm0(GetVecOfConstPtrs(chiIncr)) * m_dt / m_dpdtFactor;
    amrex::Print() << "      Before SDC " << a_sdcIter
                   << ": max relative P mismatch is " << max_corr << "\n";
  }
}

void
PeleLM::macProject(
  const TimeStamp a_time,
  const std::unique_ptr<AdvanceAdvData>& advData,
  const amrex::Vector<amrex::MultiFab*>& a_divu)
{
  BL_PROFILE("PeleLMeX::macProject()");

  const int has_divu = static_cast<int>(!a_divu.empty());

  // Get face rho inv
  auto bcRec = fetchBCRecArray(DENSITY, 1);
  amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> rho_inv(
    finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (m_incompressible != 0) {
      amrex::Real rhoInv = m_dt / (2.0 * m_rho);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        rho_inv[lev][idim].define(
          amrex::convert(grids[lev], amrex::IntVect::TheDimensionVector(idim)),
          dmap[lev], 1, 0, amrex::MFInfo(), Factory(lev));
        rho_inv[lev][idim].setVal(rhoInv);
      }
    } else {
      auto* ldata_p = getLevelDataPtr(lev, a_time);
      constexpr int doZeroVisc = 0;
      rho_inv[lev] =
        getDiffusivity(lev, DENSITY, 1, doZeroVisc, {bcRec}, ldata_p->state);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        rho_inv[lev][idim].invert(m_dt / 2.0, 0);
        rho_inv[lev][idim].FillBoundary(geom[lev].periodicity());
      }
    }
  }

  // For closed chamber, compute change in chamber pressure
  amrex::Real Sbar = 0.0;
  if ((m_closed_chamber != 0) && (m_incompressible == 0)) {
    Sbar = adjustPandDivU(advData);
  }

  if (macproj->needInitialization()) {
    amrex::LPInfo lpInfo;
    lpInfo.setMaxCoarseningLevel(m_mac_mg_max_coarsening_level);
    macproj->initProjector(lpInfo, GetVecOfArrOfConstPtrs(rho_inv));
    macproj->setDomainBC(
      getMACProjectionBC(amrex::Orientation::low),
      getMACProjectionBC(amrex::Orientation::high));
#ifdef AMREX_USE_HYPRE
    if (
      macproj->getMLMG().getBottomSolver() ==
      amrex::MLMG::BottomSolver::hypre) {
      macproj->getMLMG().setHypreOptionsNamespace(m_hypre_namespace_mac);
    }
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

#ifdef AMREX_USE_EB
  if (m_useEBinflow != 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      macproj->setEBInflowVelocity(
        lev, *getEBState(lev, VELX, AMREX_SPACEDIM, a_time));
    }
  }
#endif

  // Project
  macproj->project(m_mac_mg_rtol, m_mac_mg_atol);

  // Restore mac_divu
  if ((m_closed_chamber != 0) && (m_incompressible == 0)) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      a_divu[lev]->plus(Sbar, 0, 1);
    }
  }

  // FillBoundary umac
  // Do coarse first
  AMREX_D_TERM(advData->umac[0][0].FillBoundary(geom[0].periodicity());
               , advData->umac[0][1].FillBoundary(geom[0].periodicity());
               , advData->umac[0][2].FillBoundary(geom[0].periodicity()));
  for (int lev = 1; lev <= finest_level; ++lev) {
    // We need to fill the MAC velocities outside the fine region so we can
    // use them in the Godunov method
    const amrex::IntVect rr =
      geom[lev].Domain().size() / geom[lev - 1].Domain().size();
    create_constrained_umac_grown(
      m_nGrowMAC, &geom[lev - 1], &geom[lev],
      GetArrOfPtrs(advData->umac[lev - 1]), GetArrOfPtrs(advData->umac[lev]),
      rr);
  }
}

void
PeleLM::create_constrained_umac_grown(
  const int a_nGrow,
  const amrex::Geometry* crse_geom,
  const amrex::Geometry* fine_geom,
  const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& u_mac_crse,
  const amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM>& u_mac_fine,
  const amrex::IntVect& crse_ratio)
{
  // Divergence preserving interp
  amrex::Interpolater* mapper = &amrex::face_divfree_interp;

  // Set BCRec for Umac
  amrex::Vector<amrex::BCRec> bcrec(1);
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (crse_geom->isPeriodic(idim)) {
      bcrec[0].setLo(idim, amrex::BCType::int_dir);
      bcrec[0].setHi(idim, amrex::BCType::int_dir);
    } else {
      bcrec[0].setLo(idim, amrex::BCType::foextrap);
      bcrec[0].setHi(idim, amrex::BCType::foextrap);
    }
  }
  amrex::Array<amrex::Vector<amrex::BCRec>, AMREX_SPACEDIM> bcrecArr = {
    AMREX_D_DECL(bcrec, bcrec, bcrec)};

  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<umacFill>> crse_bndry_func(
    *crse_geom, bcrec, umacFill{});
  amrex::Array<
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<umacFill>>, AMREX_SPACEDIM>
    cbndyFuncArr = {
      AMREX_D_DECL(crse_bndry_func, crse_bndry_func, crse_bndry_func)};

  amrex::PhysBCFunct<amrex::GpuBndryFuncFab<umacFill>> fine_bndry_func(
    *fine_geom, bcrec, umacFill{});
  amrex::Array<
    amrex::PhysBCFunct<amrex::GpuBndryFuncFab<umacFill>>, AMREX_SPACEDIM>
    fbndyFuncArr = {
      AMREX_D_DECL(fine_bndry_func, fine_bndry_func, fine_bndry_func)};

  // Use piecewise constant interpolation in time, so create dummy variable for
  // time
  constexpr amrex::Real dummy = 0.;
  FillPatchTwoLevels(
    u_mac_fine, amrex::IntVect(a_nGrow), dummy, {u_mac_crse}, {dummy},
    {u_mac_fine}, {dummy}, 0, 0, 1, *crse_geom, *fine_geom, cbndyFuncArr, 0,
    fbndyFuncArr, 0, crse_ratio, mapper, bcrecArr, 0);
}

amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM>
PeleLM::getMACProjectionBC(amrex::Orientation::Side a_side)
{

  amrex::Array<amrex::LinOpBCType, AMREX_SPACEDIM> r;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (Geom(0).isPeriodic(idim)) {
      r[idim] = amrex::LinOpBCType::Periodic;
    } else {
      auto physbc = (a_side == amrex::Orientation::low) ? m_phys_bc.lo(idim)
                                                        : m_phys_bc.hi(idim);
      if (physbc == amrex::PhysBCType::outflow) {
        r[idim] = amrex::LinOpBCType::Dirichlet;
      } else {
        r[idim] = amrex::LinOpBCType::Neumann;
      }
    }
  }
  return r;
}
