#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <PeleLMeX_EF_K.H>
#include <MLGMRES.H>
#include <PeleLMeX_DiffusionOp.H>
#include <PeleLMeX_Utils.H>
#include <hydro_utils.H>

using namespace amrex;

PrecondOp*
PeleLM::getPrecondOp()
{
  if (!m_precond_op)
    m_precond_op.reset(new PrecondOp(this));
  return m_precond_op.get();
}

void
PeleLM::implicitNonLinearSolve(
  int sdcIter,
  const Real& a_dt,
  std::unique_ptr<AdvanceDiffData>& diffData,
  std::unique_ptr<AdvanceAdvData>& advData)
{
  BL_PROFILE_VAR("PeleLMeX::implicitNonLinearSolve()", implicitNonLinearSolve);

  const Real strt_time = ParallelDescriptor::second();

  //------------------------------------------------------------------------
  // Pre-solve

  // Get cell-centered charged species transport coefficients
  calcEFTransport(AmrNewTime);

  // Substepping of non-linear solve
  dtsub = a_dt / ef_substep;

  // Pass t^{n} nE/PhiV from leveldata to leveldatanlsolve
  // t^{n} have been fillpatched already
  for (int lev = 0; lev <= finest_level; ++lev) {
    // Get t^{n} data pointer
    auto ldata_p = getLevelDataPtr(lev, AmrOldTime);
    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
    MultiFab::Copy(ldataNLs_p->nlState, ldata_p->state, NE, 0, 1, m_nGrowState);
    MultiFab::Copy(
      ldataNLs_p->nlState, ldata_p->state, PHIV, 1, 1, m_nGrowState);
  }

  // Gradient of PhiV at t^{n}
  int do_avgDown = 0; // TODO or should I ?
  auto bcRecPhiV = fetchBCRecArray(PHIV, 1);
  getDiffusionOp()->computeGradient(
    getNLgradPhiVVect(), {}, // don't need the laplacian out
    GetVecOfConstPtrs(getPhiVVect(AmrOldTime)), bcRecPhiV[0], do_avgDown);

  // Stash away a copy of umac
  for (int lev = 0; lev <= finest_level; ++lev) {
    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      MultiFab::Copy(
        ldataNLs_p->umac[idim], advData->umac[lev][idim], 0, 0, 1, 0);
    }
  }

  // Setup MLGMRES
  MLGMRESSolver gmres;
  int GMRES_tot_count = 0;
  if (!m_ef_use_PETSC_direct) {
    gmres.define(this, 2, 1);
    MLJtimesVFunc jtv = &PeleLM::jTimesV;
    gmres.setJtimesV(jtv);
    MLNormFunc normF = &PeleLM::nlSolveNorm;
    gmres.setNorm(normF);
    MLPrecondFunc prec = &PeleLM::applyPrecond;
    gmres.setPrecond(prec);
  }

  //------------------------------------------------------------------------
  // Outer subcycling loop
  if (ef_substep > 1)
    Abort("Non-linear solve sub-stepping not re-implemented yet");
  int NK_tot_count = 0;
  for (int sstep = 0; sstep < ef_substep; sstep++) {

    curtime = getTime(0, AmrOldTime) + (sstep + 1) * dtsub;

    // -----------------
    // Pre-Newton
    // Set up the NL state scaling
    getNLStateScaling(nE_scale, phiV_scale);
    nE_scale = (nE_scale > 1.0e-12) ? nE_scale : 1.0;
    phiV_scale = (phiV_scale > 1.0e-12) ? phiV_scale : 1.0;
    if (ef_verbose) {
      amrex::Print() << "(" << sstep << ") ne scaling: " << nE_scale << "\n";
      amrex::Print() << "(" << sstep << ") phiV scaling: " << phiV_scale
                     << "\n";
    }
    scaleNLState(nE_scale, phiV_scale);

    // Compute the background charge distribution at curtime
    computeBGcharge(curtime, diffData, advData);
    // WriteDebugPlotFile(GetVecOfConstPtrs(getNLBGChargeVect()),"NLBGCharge");

    // Newton initial guess
    // TODO newton initial guess
    nlSolveNorm(getNLstateVect(), nl_stateNorm);

    // Initial NL residual: update residual scaling and preconditioner
    int update_scaling = 1;
    int update_precond = 1;
    nonLinearResidual(
      dtsub, getNLstateVect(), getNLresidVect(), update_scaling,
      update_precond);
    nlSolveNorm(getNLresidVect(), nl_residNorm);
    // WriteDebugPlotFile(GetVecOfConstPtrs(getNLresidVect()),"NLResInit");

    // Check for convergence
    Real scaledResnE, scaledResphiV;
    getNLResidScaling(scaledResnE, scaledResphiV);
    Real max_nlres = std::max(scaledResnE, scaledResphiV);
    if (max_nlres <= m_ef_newtonTol) {
      if (ef_verbose) {
        amrex::Print() << " No Newton iteration needed, exiting. \n";
      }
      return;
    }

    // -----------------
    // Newton iteration
    int exit_newton = 0;
    int NK_ite = 0;
    do {
      NK_ite += 1;

      // Verbose
      if (ef_verbose) {
        amrex::Print() << " Newton it: " << NK_ite << " L2**2 residual: "
                       << 0.5 * nl_residNorm * nl_residNorm
                       << ". Linf residual: " << max_nlres << "\n";
      }

      // Solve for Newton direction
      Vector<MultiFab> newtonDir(finest_level + 1);
      for (int lev = 0; lev <= finest_level; ++lev) {
        newtonDir[lev].define(
          grids[lev], dmap[lev], 2, 1, MFInfo(), Factory(lev));
        newtonDir[lev].setVal(0.0, 0, 2, 1);
      }
      if (!m_ef_use_PETSC_direct) {
        const Real S_tol = m_ef_GMRES_reltol;
        const Real S_tol_abs = m_ef_GMRES_abstol;
        GMRES_tot_count += gmres.solve(
          GetVecOfPtrs(newtonDir), getNLresidVect(), S_tol_abs, S_tol);
      } else {
      }
      // WriteDebugPlotFile(GetVecOfConstPtrs(newtonDir),"newtonDir_"+std::to_string(NK_ite));
      Real newtonDir_Norm;
      nlSolveNorm(GetVecOfPtrs(newtonDir), newtonDir_Norm);
      // Print() << " newtonDir_Norm " << newtonDir_Norm << "\n";

      // Linesearch & update state TODO
      updateNLState(GetVecOfPtrs(newtonDir));
      nlSolveNorm(getNLstateVect(), nl_stateNorm);
      update_scaling = 0;
      update_precond = 1;
      nonLinearResidual(
        dtsub, getNLstateVect(), getNLresidVect(), update_scaling,
        update_precond);
      nlSolveNorm(getNLresidVect(), nl_residNorm);
      getNLResidScaling(scaledResnE, scaledResphiV);
      max_nlres = std::max(scaledResnE, scaledResphiV);
      // WriteDebugPlotFile(GetVecOfConstPtrs(getNLstateVect()),"NLState_"+std::to_string(NK_ite));
      // WriteDebugPlotFile(GetVecOfConstPtrs(getNLresidVect()),"NLResid_"+std::to_string(NK_ite));

      // Exit condition
      exit_newton = testExitNewton(NK_ite, max_nlres, newtonDir_Norm);
    } while (!exit_newton);
    NK_tot_count += NK_ite;

    // -----------------
    // Post-Newton
    // Increment the forcing term
    incrementElectronForcing(sstep, advData);
    // Unscale nl_state and if not last subcycle update 'old' state
    for (int lev = 0; lev <= finest_level; ++lev) {
      m_leveldatanlsolve[lev]->nlState.mult(nE_scale, 0, 1, 0);
      m_leveldatanlsolve[lev]->nlState.mult(phiV_scale, 1, 1, 0);
    }
  }

  // Update the state
  for (int lev = 0; lev <= finest_level; ++lev) {
    // Get t^{n} data pointer
    auto ldata_p = getLevelDataPtr(lev, AmrNewTime);
    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
    int nGrowNL = 0;
    MultiFab::Copy(ldata_p->state, ldataNLs_p->nlState, 0, NE, 1, nGrowNL);
    MultiFab::Copy(ldata_p->state, ldataNLs_p->nlState, 1, PHIV, 1, nGrowNL);
  }

  if (ef_verbose) {
    Real run_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(
      run_time, ParallelDescriptor::IOProcessorNumber());
    if (!m_ef_use_PETSC_direct) {
      Real avgGMRES = (float)GMRES_tot_count / (float)NK_tot_count;
      amrex::Print() << "  [" << sdcIter << "] dt: " << a_dt
                     << " - Avg GMRES/Newton: " << avgGMRES << "\n";
    }
    amrex::Print() << "  >> PeleLMeX::implicitNLSolve() " << run_time << "\n";
  }

  // VisMF::Write(advData->Forcing[0],"ForcingNE");
}

int
PeleLM::testExitNewton(
  int newtonIter, const Real& max_res, const Real& norm_NewtonDir)
{
  int exit = 0;
  if (max_res <= m_ef_newtonTol || norm_NewtonDir <= 1e-11) {
    exit = 1;
    if (ef_verbose) {
      amrex::Print() << " Newton iterations converged: \n";
      amrex::Print() << " Final Newton L2**2 res norm : "
                     << 0.5 * nl_residNorm * nl_residNorm << "\n";
      amrex::Print() << " Final Newton Linf res norm : " << max_res << "\n";
    }
  }

  if (newtonIter >= m_ef_maxNewtonIter && exit == 0) {
    exit = 1;
    amrex::Print()
      << " WARNING: Max Newton iteration reached without convergence !!! \n";
    amrex::Print() << " Final Newton L2**2 res norm : "
                   << 0.5 * nl_residNorm * nl_residNorm << "\n";
    amrex::Print() << " Final Newton Linf res norm : " << max_res << "\n";
  }
  return exit;
}

void
PeleLM::updateNLState(const Vector<MultiFab*>& a_update)
{
  // AverageDown the newton direction
  /*
  for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
     EB_average_down(*a_update[lev],
                     *a_update[lev-1],
                     0,2,refRatio(lev-1));
#else
     average_down(*a_update[lev],
                  *a_update[lev-1],
                  0,2,refRatio(lev-1));
#endif
  }
  */
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev); // NL data
    ldataNLs_p->nlState.plus(*a_update[lev], 0, 2, 0);
  }

  // AverageDown the new state
  for (int lev = finest_level; lev > 0; --lev) {
    auto ldataNLsFine_p = getLevelDataNLSolvePtr(lev);     // NL data
    auto ldataNLsCrse_p = getLevelDataNLSolvePtr(lev - 1); // NL data
#ifdef AMREX_USE_EB
    EB_average_down(
      ldataNLsFine_p->nlState, ldataNLsCrse_p->nlState, 0, 2,
      refRatio(lev - 1));
#else
    average_down(
      ldataNLsFine_p->nlState, ldataNLsCrse_p->nlState, 0, 2,
      refRatio(lev - 1));
#endif
  }

  // FillBoundary NLState
  /*
  for (int lev = 0; lev <= finest_level; ++lev) {
     auto ldataNLs_p = getLevelDataNLSolvePtr(lev);   // NL data
     ldataNLs_p->nlState.FillBoundary(geom[lev].periodicity());
  }
  */

  // Need to fillpatch the NL state
  // First unscale NLState
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev); // NL data
    ldataNLs_p->nlState.mult(nE_scale, 0, 1, m_nGrowState);
    ldataNLs_p->nlState.mult(phiV_scale, 1, 1, m_nGrowState);
  }

  // FillPatch
  Vector<MultiFab> nEState;
  Vector<MultiFab> phiVState;
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev); // NL data
    nEState.emplace_back(ldataNLs_p->nlState, amrex::make_alias, 0, 1);
    phiVState.emplace_back(ldataNLs_p->nlState, amrex::make_alias, 1, 1);
  }
  fillPatchNLnE(m_cur_time, GetVecOfPtrs(nEState), m_nGrowState);
  fillPatchNLphiV(m_cur_time, GetVecOfPtrs(phiVState), m_nGrowState);

  // Rescale
  scaleNLState(nE_scale, phiV_scale);
}

void
PeleLM::incrementElectronForcing(
  int a_sstep, std::unique_ptr<AdvanceAdvData>& advData)
{
  for (int lev = 0; lev <= finest_level; ++lev) {

    auto ldata_p = getLevelDataPtr(lev, AmrOldTime); // Old time electron
    auto ldataR_p = getLevelDataReactPtr(lev);       // Reaction
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);   // NL data

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldata_p->state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& nE_o = ldata_p->state.const_array(mfi, NE);
      auto const& nE_n = ldataNLs_p->nlState.const_array(mfi);
      auto const& I_R_nE = ldataR_p->I_R.const_array(mfi, NUM_SPECIES);
      auto const& FnE = advData->Forcing[lev].array(mfi, NUM_SPECIES + 1);
      Real scaling = nE_scale;
      Real dtinv = 1.0 / dtsub;
      amrex::ParallelFor(
        bx, [nE_o, nE_n, I_R_nE, FnE, dtinv, scaling,
             a_sstep] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (a_sstep == 0) {
            FnE(i, j, k) = (nE_n(i, j, k) * scaling - nE_o(i, j, k)) * dtinv -
                           I_R_nE(i, j, k);
          } else {
            FnE(i, j, k) += (nE_n(i, j, k) * scaling - nE_o(i, j, k)) * dtinv -
                            I_R_nE(i, j, k);
          }
        });
    }
  }
}

void
PeleLM::computeBGcharge(
  const Real& a_time,
  std::unique_ptr<AdvanceDiffData>& diffData,
  std::unique_ptr<AdvanceAdvData>& advData)
{
  // Get integration dt
  Real dt_int = a_time - getTime(0, AmrOldTime);

  for (int lev = 0; lev <= finest_level; ++lev) {

    // Get data pointers
    auto ldata_p = getLevelDataPtr(lev, AmrOldTime); // Old time species
    auto ldataR_p = getLevelDataReactPtr(lev);       // Reaction
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);   // NL data

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldataNLs_p->backgroundCharge, TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& rhoYold = ldata_p->state.const_array(mfi, FIRSTSPEC);
      auto const& adv_arr = advData->AofS[lev].const_array(mfi, FIRSTSPEC);
      auto const& dn_arr = diffData->Dn[lev].const_array(mfi);
      auto const& dnp1_arr = diffData->Dnp1[lev].const_array(mfi);
      auto const& dhat_arr = diffData->Dhat[lev].const_array(mfi);
      auto const& rhoYdot = ldataR_p->I_R.const_array(mfi);
      auto const& charge = ldataNLs_p->backgroundCharge.array(mfi);
      Real factor = 1.0 / elemCharge;
      amrex::ParallelFor(
        bx,
        [dt_int, rhoYold, adv_arr, dn_arr, dnp1_arr, dhat_arr, rhoYdot, charge,
         factor, zk = zk] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          charge(i, j, k) = 0.0;
          for (int n = 0; n < NUM_SPECIES; n++) {
            Real rhoYprov =
              rhoYold(i, j, k, n) +
              dt_int * (adv_arr(i, j, k, n) +
                        0.5 * (dn_arr(i, j, k, n) - dnp1_arr(i, j, k, n)) +
                        dhat_arr(i, j, k, n) + rhoYdot(i, j, k, n));
            rhoYprov = amrex::max(rhoYprov, 0.0);
            charge(i, j, k) += zk[n] * rhoYprov;
          }
          charge(i, j, k) *= factor;
        });
    }
  }
}

void
PeleLM::nonLinearResidual(
  const Real& a_dt,
  const Vector<MultiFab*>& a_nlstate,
  const Vector<MultiFab*>& a_nlresid,
  int updateScaling,
  int updatePrecond)
{
  // Get unscaled copy of the NL state
  Vector<MultiFab> nE(finest_level + 1);
  Vector<MultiFab> phiV(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    nE[lev].define(
      grids[lev], dmap[lev], 1, m_nGrowState, MFInfo(), Factory(lev));
    phiV[lev].define(
      grids[lev], dmap[lev], 1, m_nGrowState, MFInfo(), Factory(lev));
    MultiFab::Copy(nE[lev], *a_nlstate[lev], 0, 0, 1, m_nGrowState);
    nE[lev].mult(nE_scale, 0, 1, m_nGrowState);
    MultiFab::Copy(phiV[lev], *a_nlstate[lev], 1, 0, 1, m_nGrowState);
    phiV[lev].mult(phiV_scale, 0, 1, m_nGrowState);
  }

  // FillPatch nE state here
  fillPatchNLnE(m_cur_time, GetVecOfPtrs(nE), m_nGrowState);

  // Get L(phiV) and Grad(phiV)
  Vector<MultiFab> laplacian(finest_level + 1);
  Vector<Array<MultiFab, AMREX_SPACEDIM>> gradPhiVCur(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    laplacian[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      const auto& fba =
        amrex::convert(grids[lev], IntVect::TheDimensionVector(idim));
      gradPhiVCur[lev][idim].define(
        fba, dmap[lev], 1, 0, MFInfo(), Factory(lev));
    }
  }
  int do_avgDown = 0; // TODO or should I ?
  auto bcRecPhiV = fetchBCRecArray(PHIV, 1);
  getDiffusionOp()->computeGradient(
    GetVecOfArrOfPtrs(gradPhiVCur), GetVecOfPtrs(laplacian),
    GetVecOfConstPtrs(phiV), bcRecPhiV[0], do_avgDown);

  // Get nE diffusion term
  Vector<MultiFab> diffnE(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    diffnE[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
  }
  auto bcRecnE = fetchBCRecArray(NE, 1);
  getDiffusionOp()->computeDiffLap(
    GetVecOfPtrs(diffnE), 0, GetVecOfConstPtrs(nE), 0,
    GetVecOfConstPtrs(getnEDiffusivityVect(AmrNewTime)), 0, bcRecnE, 1);
  // WriteDebugPlotFile(GetVecOfConstPtrs(diffnE),"diffnE");
  // VisMF::Write(nE[0],"nEForDiffnlResid");
  // VisMF::Write(diffnE[0],"diffnEnlResid");

  // Get nE advection term
  Vector<MultiFab> advnE(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    advnE[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), Factory(lev));
  }
  getAdvectionTerm(
    GetVecOfConstPtrs(nE), GetVecOfPtrs(advnE),
    GetVecOfArrOfConstPtrs(gradPhiVCur));
  // WriteDebugPlotFile(GetVecOfConstPtrs(advnE),"advnE");

  // Assemble non-linear residual
  // res(ne(:)) = dt * ( diff(:) + conv(:) + I_R(:) ) - ( ne(:) - ne_old(:) )
  // res(phiv(:)) = \Sum z_k * \tilde Y_k / q_e - ne + Lapl_PhiV
  for (int lev = 0; lev <= finest_level; ++lev) {

    auto ldataOld_p = getLevelDataPtr(lev, AmrOldTime); // Old time
    auto ldataR_p = getLevelDataReactPtr(lev);          // Reaction

    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);

    // Init the ghostcells too
    a_nlresid[lev]->setVal(0.0);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldataNLs_p->nlResid, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& I_R_nE = ldataR_p->I_R.const_array(mfi, NUM_SPECIES);
      auto const& lapPhiV = laplacian[lev].const_array(mfi);
      auto const& ne_diff = diffnE[lev].const_array(mfi);
      auto const& ne_adv = advnE[lev].const_array(mfi);
      auto const& ne_curr = nE[lev].const_array(mfi);
      auto const& ne_old = ldataOld_p->state.const_array(mfi, NE);
      auto const& charge = ldataNLs_p->backgroundCharge.const_array(mfi);
      auto const& res_nE = a_nlresid[lev]->array(mfi, 0);
      auto const& res_phiV = a_nlresid[lev]->array(mfi, 1);
      Real scalLap = eps0 * epsr / elemCharge;
      amrex::ParallelFor(
        bx, [ne_curr, ne_old, lapPhiV, I_R_nE, ne_diff, ne_adv, charge, res_nE,
             res_phiV, a_dt,
             scalLap] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          res_nE(i, j, k) =
            ne_old(i, j, k) - ne_curr(i, j, k) +
            a_dt * (ne_diff(i, j, k) + ne_adv(i, j, k) + I_R_nE(i, j, k));
          res_phiV(i, j, k) =
            lapPhiV(i, j, k) * scalLap - ne_curr(i, j, k) + charge(i, j, k);
          res_nE(i, j, k) *= -1.0;   // NLresidual is -RHS
          res_phiV(i, j, k) *= -1.0; // NLresidual is -RHS
        });
    }
  }
  // WriteDebugPlotFile(GetVecOfConstPtrs(a_nlresid),"UnscalednlResid");

  //   /*
  // AverageDown the NLresiduals
  for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
    EB_average_down(
      *a_nlresid[lev], *a_nlresid[lev - 1], 0, 2, refRatio(lev - 1));
#else
    average_down(*a_nlresid[lev], *a_nlresid[lev - 1], 0, 2, refRatio(lev - 1));
#endif
  }
  //   */

  // Update the residual scale if necessary
  if (updateScaling) {
    getNLResidScaling(FnE_scale, FphiV_scale);
    FnE_scale = (FnE_scale > 1.0e-12) ? FnE_scale : 1.0;
    FphiV_scale = (FphiV_scale > 1.0e-12) ? FphiV_scale : 1.0;
    if (ef_verbose) {
      amrex::Print() << "    F(ne) scaling: " << FnE_scale << "\n";
      amrex::Print() << "    F(PhiV) scaling: " << FphiV_scale << "\n";
    }
  }

  // Apply scaling
  scaleNLResid(a_nlresid, FnE_scale, FphiV_scale);

  // Update the preconditioning LinOps
  if (updatePrecond && !m_ef_use_PETSC_direct) {
    setUpPrecond(a_dt, GetVecOfConstPtrs(nE));
  }
}

void
PeleLM::getAdvectionTerm(
  const Vector<const MultiFab*>& a_nE,
  const Vector<MultiFab*>& a_advTerm,
  const Vector<Array<const MultiFab*, AMREX_SPACEDIM>>& a_gPhiVCur)
{

  // Get advection fluxes on all levels
  int nGrow = 0;
  Vector<Array<MultiFab, AMREX_SPACEDIM>> fluxes(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    const auto& ba = grids[lev];
    const auto& factory = Factory(lev);
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      fluxes[lev][idim].define(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
        nGrow, MFInfo(), factory);
    }
  }

  // nE BCRec
  auto bcRecnE = fetchBCRecArray(NE, 1);

  for (int lev = 0; lev <= finest_level; ++lev) {

    // Get t^{n} data pointer
    auto ldata_p = getLevelDataPtr(lev, AmrNewTime);

    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);

    // Get the face centered electron mobility
    int doZeroVisc = 0;
    Array<MultiFab, AMREX_SPACEDIM> mobE_ec =
      getDiffusivity(lev, 0, 1, doZeroVisc, bcRecnE, ldata_p->mobE_cc);

    // Get the electron effective velocity
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNLs_p->uEffnE[idim], TilingIfNotGPU());
           mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto const& ueff = ldataNLs_p->uEffnE[idim].array(mfi);
        auto const& umac = ldataNLs_p->umac[idim].const_array(mfi);
        auto const& gphi_c = a_gPhiVCur[lev][idim]->const_array(mfi);
        auto const& gphi_o = ldataNLs_p->gPhiVOld[idim].const_array(mfi);
        auto const& kappa_e = mobE_ec[idim].const_array(mfi);
        amrex::ParallelFor(
          bx, [umac, ueff, gphi_c, gphi_o,
               kappa_e] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            ueff(i, j, k) =
              umac(i, j, k) - kappa_e(i, j, k) * -1.0 * 0.5 *
                                (gphi_c(i, j, k) + gphi_o(i, j, k));
          });
      }
    }
  }

  // Average down the Ueff
  /*
  for (int lev = finest_level; lev > 0; --lev) {
     auto ldataNLFine_p = getLevelDataNLSolvePtr(lev);
     auto ldataNLCrse_p = getLevelDataNLSolvePtr(lev-1);
#ifdef AMREX_USE_EB
     EB_average_down_faces(GetArrOfConstPtrs(ldataNLFine_p->uEffnE),
                           GetArrOfPtrs(ldataNLCrse_p->uEffnE),
                           refRatio(lev-1),nGrow);
#else
     average_down_faces(GetArrOfConstPtrs(ldataNLFine_p->uEffnE),
                        GetArrOfPtrs(ldataNLCrse_p->uEffnE),
                        refRatio(lev-1),nGrow);
#endif
  }
  */

  for (int lev = 0; lev <= finest_level; ++lev) {
    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
    // Compute advective fluxes
    getAdvectionFluxes(
      lev, GetArrOfPtrs(fluxes[lev]), *a_nE[lev],
      GetArrOfConstPtrs(ldataNLs_p->uEffnE), bcRecnE[0]);
    // getAdvectionFluxesMOL(lev, GetArrOfPtrs(fluxes[lev]), *a_nE[lev],
    // GetArrOfConstPtrs(ldataNLs_p->uEffnE), bcRecnE[0]);
    // VisMF::Write(fluxes[lev][1],"FluxYNewbefAvgDown_lev"+std::to_string(lev));
  }

  // Average down the fluxes
  for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
    EB_average_down_faces(
      GetArrOfConstPtrs(fluxes[lev]), GetArrOfPtrs(fluxes[lev - 1]),
      refRatio(lev - 1), geom[lev - 1]);
#else
    average_down_faces(
      GetArrOfConstPtrs(fluxes[lev]), GetArrOfPtrs(fluxes[lev - 1]),
      refRatio(lev - 1), geom[lev - 1]);
#endif
  }

  // Compute divergence
  int intensiveFluxes = 1;
  fluxDivergence(
    a_advTerm, 0, GetVecOfArrOfPtrs(fluxes), 0, 1, intensiveFluxes, -1.0);
}

void
PeleLM::getAdvectionFluxesMOL(
  int lev,
  const Array<MultiFab*, AMREX_SPACEDIM>& a_fluxes,
  const MultiFab& a_nE,
  const Array<const MultiFab*, AMREX_SPACEDIM>& a_ueff,
  BCRec bcrec)
{
  auto bcRec_d = convertToDeviceVector(Vector<BCRec>(1, bcrec));
  auto AdvType = fetchAdvTypeArray(FIRSTSPEC, 1);
  auto AdvType_d = convertToDeviceVector(AdvType);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox edgstate[AMREX_SPACEDIM];
    for (MFIter mfi(a_nE, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

      Box const& bx = mfi.tilebox();
      AMREX_D_TERM(const Box& xbx = mfi.grownnodaltilebox(0, 0);
                   , const Box& ybx = mfi.grownnodaltilebox(1, 0);
                   , const Box& zbx = mfi.grownnodaltilebox(2, 0));

      AMREX_D_TERM(auto const& fx = a_fluxes[0]->array(mfi, 0);
                   , auto const& fy = a_fluxes[1]->array(mfi, 0);
                   , auto const& fz = a_fluxes[2]->array(mfi, 0);)
      AMREX_D_TERM(auto const& ueff = a_ueff[0]->const_array(mfi);
                   , auto const& veff = a_ueff[1]->const_array(mfi);
                   , auto const& weff = a_ueff[2]->const_array(mfi);)
      AMREX_D_TERM(edgstate[0].resize(xbx, 1);, edgstate[1].resize(ybx, 1);
                   , edgstate[2].resize(zbx, 1));
      AMREX_D_TERM(Array4<Real> xstate = edgstate[0].array();
                   , Array4<Real> ystate = edgstate[1].array();
                   , Array4<Real> zstate = edgstate[2].array());
      AMREX_D_TERM(Elixir eli_edgex = edgstate[0].elixir();
                   , Elixir eli_edgey = edgstate[1].elixir();
                   , Elixir eli_edgez = edgstate[2].elixir());

      auto const& nE_arr = a_nE.const_array(mfi);
      auto const& divu_arr = a_nE.const_array(mfi);
      auto const& force_arr = a_nE.const_array(mfi);

      bool is_velocity = false;
      bool fluxes_are_area_weighted = false;
      bool knownEdgeState = false;
      std::string mol = "MOL";
      HydroUtils::ComputeFluxesOnBoxFromState(
        bx, 1, mfi, nE_arr, AMREX_D_DECL(fx, fy, fz),
        AMREX_D_DECL(xstate, ystate, zstate), knownEdgeState,
        AMREX_D_DECL(ueff, veff, weff), divu_arr, force_arr, geom[lev], m_dt,
        {bcrec}, bcRec_d.dataPtr(), AdvType_d.dataPtr(),
#ifdef AMREX_USE_EB
        EBFactory(lev),
#endif
        0, 0, is_velocity, fluxes_are_area_weighted, mol);
    }
  }
}

void
PeleLM::getAdvectionFluxes(
  int lev,
  const Array<MultiFab*, AMREX_SPACEDIM>& a_fluxes,
  const MultiFab& a_nE,
  const Array<const MultiFab*, AMREX_SPACEDIM>& a_ueff,
  BCRec bcrec)
{
  const Box& domain = geom[lev].Domain();

  int order = m_nEAdvOrder;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  {
    FArrayBox edgstate[AMREX_SPACEDIM];
    for (MFIter mfi(a_nE, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Box const& bx = mfi.tilebox();
      AMREX_D_TERM(const Box& xbx = surroundingNodes(bx, 0);
                   , const Box& ybx = surroundingNodes(bx, 1);
                   , const Box& zbx = surroundingNodes(bx, 2));

      // data arrays
      auto const& ne_arr = a_nE.const_array(mfi);
      AMREX_D_TERM(Array4<Real> xflux = a_fluxes[0]->array(mfi);
                   , Array4<Real> yflux = a_fluxes[1]->array(mfi);
                   , Array4<Real> zflux = a_fluxes[2]->array(mfi));
      AMREX_D_TERM(Array4<Real const> u = a_ueff[0]->const_array(mfi);
                   , Array4<Real const> v = a_ueff[1]->const_array(mfi);
                   , Array4<Real const> w = a_ueff[2]->const_array(mfi););
      AMREX_D_TERM(edgstate[0].resize(xbx, 1);, edgstate[1].resize(ybx, 1);
                   , edgstate[2].resize(zbx, 1));
      AMREX_D_TERM(Array4<Real> xstate = edgstate[0].array();
                   , Array4<Real> ystate = edgstate[1].array();
                   , Array4<Real> zstate = edgstate[2].array());
      AMREX_D_TERM(Elixir xstate_eli = edgstate[0].elixir();
                   , Elixir ystate_eli = edgstate[1].elixir();
                   , Elixir zstate_eli = edgstate[2].elixir());

      // Predict edge states
      // X
      {
        // BCs
        const Box& edomain = surroundingNodes(domain, 0);
        const auto bc_lo = bcrec.lo(0);
        const auto bc_hi = bcrec.hi(0);

        amrex::ParallelFor(
          xbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int idx[3] = {i, j, k};
            bool on_lo =
              ((bc_lo == BCType::ext_dir) && (idx[0] <= edomain.smallEnd(0)));
            bool on_hi =
              ((bc_hi == BCType::ext_dir) && (idx[0] >= edomain.bigEnd(0)));
            if (order == 1) {
              xstate(i, j, k) =
                ef_edge_state_extdir(i, j, k, 0, on_lo, on_hi, ne_arr, u);
            } else if (order == 2) {
              bool extdir_or_ho_lo =
                (bc_lo == BCType::ext_dir) || (bc_lo == BCType::hoextrap);
              bool extdir_or_ho_hi =
                (bc_hi == BCType::ext_dir) || (bc_hi == BCType::hoextrap);
              xstate(i, j, k) = ef_edge_state_2ndO_extdir(
                i, j, k, 0, on_lo, on_hi, extdir_or_ho_lo, extdir_or_ho_hi,
                domain.smallEnd(0), domain.bigEnd(0), ne_arr, u, xbx);
            }
          });
      }
#if (AMREX_SPACEDIM > 1)
      // Y
      {
        // BCs
        const Box& edomain = surroundingNodes(domain, 1);
        const auto bc_lo = bcrec.lo(1);
        const auto bc_hi = bcrec.hi(1);

        amrex::ParallelFor(
          ybx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int idx[3] = {i, j, k};
            bool on_lo =
              ((bc_lo == BCType::ext_dir) && (idx[1] <= edomain.smallEnd(1)));
            bool on_hi =
              ((bc_hi == BCType::ext_dir) && (idx[1] >= edomain.bigEnd(1)));
            if (order == 1) {
              ystate(i, j, k) =
                ef_edge_state_extdir(i, j, k, 1, on_lo, on_hi, ne_arr, v);
            } else if (order == 2) {
              bool extdir_or_ho_lo =
                (bc_lo == BCType::ext_dir) || (bc_lo == BCType::hoextrap);
              bool extdir_or_ho_hi =
                (bc_hi == BCType::ext_dir) || (bc_hi == BCType::hoextrap);
              ystate(i, j, k) = ef_edge_state_2ndO_extdir(
                i, j, k, 1, on_lo, on_hi, extdir_or_ho_lo, extdir_or_ho_hi,
                domain.smallEnd(1), domain.bigEnd(1), ne_arr, v, ybx);
            }
          });
      }
#if (AMREX_SPACEDIM == 3)
      // Z
      {
        // BCs
        const Box& edomain = surroundingNodes(domain, 2);
        const auto bc_lo = bcrec.lo(2);
        const auto bc_hi = bcrec.hi(2);

        amrex::ParallelFor(
          zbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            int idx[3] = {i, j, k};
            bool on_lo =
              ((bc_lo == BCType::ext_dir) && (idx[2] <= edomain.smallEnd(2)));
            bool on_hi =
              ((bc_hi == BCType::ext_dir) && (idx[2] >= edomain.bigEnd(2)));
            if (order == 1) {
              zstate(i, j, k) =
                ef_edge_state_extdir(i, j, k, 2, on_lo, on_hi, ne_arr, w);
            } else if (order == 2) {
              bool extdir_or_ho_lo =
                (bc_lo == BCType::ext_dir) || (bc_lo == BCType::hoextrap);
              bool extdir_or_ho_hi =
                (bc_hi == BCType::ext_dir) || (bc_hi == BCType::hoextrap);
              zstate(i, j, k) = ef_edge_state_2ndO_extdir(
                i, j, k, 2, on_lo, on_hi, extdir_or_ho_lo, extdir_or_ho_hi,
                domain.smallEnd(2), domain.bigEnd(2), ne_arr, w, zbx);
            }
          });
      }
#endif
#endif

      // Computing fluxes
      amrex::ParallelFor(
        xbx, [u, xstate, xflux] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          xflux(i, j, k) = u(i, j, k) * xstate(i, j, k);
        });
#if (AMREX_SPACEDIM > 1)
      amrex::ParallelFor(
        ybx, [v, ystate, yflux] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          yflux(i, j, k) = v(i, j, k) * ystate(i, j, k);
        });
#if (AMREX_SPACEDIM == 3)
      amrex::ParallelFor(
        zbx, [w, zstate, zflux] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          zflux(i, j, k) = w(i, j, k) * zstate(i, j, k);
        });
#endif
#endif
    }
  }
}

void
PeleLM::setUpPrecond(const Real& a_dt, const Vector<const MultiFab*>& a_nE)
{
  BL_PROFILE("PeleLMeX::setUpPrecond()");

  // Update LinearOps defs if needed -> done internally by the getPrecondOp()
  // func

  // nE BCRec
  auto bcRecnE = fetchBCRecArray(NE, 1);
  auto bcRecPhiV = fetchBCRecArray(PHIV, 1);

  //--------------------------------------------------------------------------
  // Set diff/drift operator
  getPrecondOp()->setDiffOpScalars(
    -nE_scale / FnE_scale, -a_dt * nE_scale / FnE_scale,
    a_dt * nE_scale / FnE_scale);
  Real omega = m_ABecCecOmega;
  getPrecondOp()->setDiffOpRelaxation(omega);
  getPrecondOp()->setDiffOpBCs(bcRecnE[0]);

  // Get and set coefficients on each level
  for (int lev = 0; lev <= finest_level; ++lev) {

    // Get level data pointer
    auto ldata_p = getLevelDataPtr(lev, AmrNewTime);

    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);

    int doZeroVisc = 0;
    Array<MultiFab, AMREX_SPACEDIM> diffE_ec =
      getDiffusivity(lev, 0, 1, doZeroVisc, bcRecnE, ldata_p->diffE_cc);

    getPrecondOp()->setDiffOpACoeff(lev, 1.0);
    getPrecondOp()->setDiffOpBCoeff(lev, GetArrOfConstPtrs(diffE_ec));
    getPrecondOp()->setDiffOpCCoeff(lev, GetArrOfConstPtrs(ldataNLs_p->uEffnE));
  }

  // Stilda approx first-order
  Vector<MultiFab> diagDiffOp(finest_level + 1);
  if (m_ef_PC_approx == 2) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      diagDiffOp[lev].define(grids[lev], dmap[lev], 1, 1);
      getPrecondOp()->getDiffOpDiagonal(lev, diagDiffOp[lev]);
      diagDiffOp[lev].mult(FnE_scale / nE_scale);
    }
    fillPatchExtrap(m_cur_time, GetVecOfPtrs(diagDiffOp), 1);
  }

  //--------------------------------------------------------------------------
  // Set Stilda and drift operators

  // Set scalars
  getPrecondOp()->setDriftOpScalars(0.0, 0.5 * phiV_scale / FnE_scale * a_dt);
  if (m_ef_PC_approx == 1 || m_ef_PC_approx == 2) {
    getPrecondOp()->setStildaOpScalars(0.0, -1.0);
  } else if (m_ef_PC_approx == 3 || m_ef_PC_approx == 4) {
    getPrecondOp()->setStildaOpScalars(0.0, 1.0);
  }

  // Set BCs
  getPrecondOp()->setDriftOpBCs(bcRecPhiV[0]);
  getPrecondOp()->setStildaOpBCs(bcRecPhiV[0]);

  // Get and set coefficients on each level
  for (int lev = 0; lev <= finest_level; ++lev) {

    // Get nl solve data pointer
    auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
    auto ldata_p = getLevelDataPtr(lev, AmrNewTime);

    // CC neKe values
    MultiFab nEKe(grids[lev], dmap[lev], 1, 1);

    // CC Schur neKe values
    MultiFab Schur_nEKe;
    if (m_ef_PC_approx == 2) {
      Schur_nEKe.define(grids[lev], dmap[lev], 1, 1);
    }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(nEKe, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& gbx = mfi.growntilebox();
      auto const& neke = nEKe.array(mfi);
      auto const& kappaE = ldata_p->mobE_cc.array(mfi);
      auto const& ne_arr = a_nE[lev]->const_array(mfi);
      auto const& Schur =
        (m_ef_PC_approx == 2) ? Schur_nEKe.array(mfi) : nEKe.array(mfi);
      auto const& diffOp_diag =
        (m_ef_PC_approx == 2) ? diagDiffOp[lev].array(mfi) : nEKe.array(mfi);
      int do_Schur = (m_ef_PC_approx == 2) ? 1 : 0;
      amrex::ParallelFor(
        gbx, [neke, ne_arr, kappaE, Schur, diffOp_diag, a_dt,
              do_Schur] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          neke(i, j, k) = kappaE(i, j, k) * ne_arr(i, j, k);
          if (do_Schur) {
            Schur(i, j, k) = -a_dt * 0.5 * neke(i, j, k) / diffOp_diag(i, j, k);
          }
        });
    }

    // Upwinded edge neKe values
    Array<MultiFab, AMREX_SPACEDIM> neKe_ec = getUpwindedEdge(
      lev, 0, 1, bcRecnE, nEKe, GetArrOfConstPtrs(ldataNLs_p->uEffnE));

    // Set drift Op coefficients
    getPrecondOp()->setDriftOpBCoeff(lev, GetArrOfConstPtrs(neKe_ec));

    // Set Stilda Op coefficients
    if (m_ef_PC_approx == 1) { // Assuming identity of the inverse of DiffOp
      // Add Stilda pieces
      Real scalLap = eps0 * epsr / elemCharge;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        neKe_ec[idim].mult(0.5 * a_dt, 0, 1);
        neKe_ec[idim].plus(scalLap, 0, 1);
      }
      getPrecondOp()->setStildaOpBCoeff(lev, GetArrOfConstPtrs(neKe_ec));
    } else if (m_ef_PC_approx == 2) { // Assuming inverse of the diag of DiffOp
      // Upwinded Schur edge neKe values
      Array<MultiFab, AMREX_SPACEDIM> Schur_neKe_ec = getUpwindedEdge(
        lev, 0, 1, bcRecnE, Schur_nEKe, GetArrOfConstPtrs(ldataNLs_p->uEffnE));
      Real scalLap = eps0 * epsr / elemCharge;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        Schur_neKe_ec[idim].plus(scalLap, 0, 1);
      }
      getPrecondOp()->setStildaOpBCoeff(lev, GetArrOfConstPtrs(Schur_neKe_ec));
    } else {
      Abort("Preconditioner option /= 1 or 2 not available yet");
    }
  }
}

Array<MultiFab, AMREX_SPACEDIM>
PeleLM::getUpwindedEdge(
  int lev,
  int edge_comp,
  int ncomp,
  Vector<BCRec> bcrec,
  const MultiFab& ccMF,
  const Array<const MultiFab*, AMREX_SPACEDIM>& ecVel)
{
  AMREX_ASSERT(bcrec.size() >= ncomp);
  AMREX_ASSERT(ccMF.nComp() >= edge_comp + ncomp);
  const auto& ba = ccMF.boxArray();
  const auto& dm = ccMF.DistributionMap();
  const auto& factory = ccMF.Factory();
  Array<MultiFab, AMREX_SPACEDIM> ecMFs{AMREX_D_DECL(
    MultiFab(
      amrex::convert(ba, IntVect::TheDimensionVector(0)), dm, ncomp, 0,
      MFInfo(), factory),
    MultiFab(
      amrex::convert(ba, IntVect::TheDimensionVector(1)), dm, ncomp, 0,
      MFInfo(), factory),
    MultiFab(
      amrex::convert(ba, IntVect::TheDimensionVector(2)), dm, ncomp, 0,
      MFInfo(), factory))};

#ifdef AMREX_USE_EB
  // TODO cen2edg_upwind EB
#else
  const Box& domain = geom[lev].Domain();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(ccMF, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      const Box ebx = mfi.nodaltilebox(idim);
      const Box& edomain = amrex::surroundingNodes(domain, idim);
      const auto& ccVal = ccMF.const_array(mfi, edge_comp);
      const auto& ecUeff = ecVel[idim]->const_array(mfi);
      const auto& ecVal = ecMFs[idim].array(mfi);
      const auto bc_lo = bcrec[0].lo(idim);
      const auto bc_hi = bcrec[0].hi(idim);
      amrex::ParallelFor(
        ebx, [idim, bc_lo, bc_hi, ccVal, ecUeff, ecVal,
              edomain] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          int idx[3] = {i, j, k};
          bool on_lo =
            ((bc_lo == amrex::BCType::ext_dir) &&
             (idx[idim] <= edomain.smallEnd(idim)));
          bool on_hi =
            ((bc_hi == amrex::BCType::ext_dir) &&
             (idx[idim] >= edomain.bigEnd(idim)));
          cen2edg_upwind(i, j, k, idim, 1, on_lo, on_hi, ecUeff, ccVal, ecVal);
        });
    }
  }
#endif

  return ecMFs;
}

void
PeleLM::jTimesV(const Vector<MultiFab*>& a_v, const Vector<MultiFab*>& a_Jv)
{
  Real vNorm;
  nlSolveNorm(a_v, vNorm);

  // v is zero, Jv is zero and return
  if (vNorm == 0.0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      a_Jv[lev]->setVal(0.0);
    }
    return;
  }

  // TODO: only one-sided difference for now
  Real delta_pert =
    m_ef_lambda_jfnk * (m_ef_lambda_jfnk + nl_stateNorm / vNorm);

  if (m_ef_diffT_jfnk == 1) {
    Vector<MultiFab> statePert(finest_level + 1);
    Vector<MultiFab> residPert(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
      statePert[lev].define(
        grids[lev], dmap[lev], 2, m_nGrowState, MFInfo(), Factory(lev));
      residPert[lev].define(
        grids[lev], dmap[lev], 2, m_nGrowState, MFInfo(), Factory(lev));
      MultiFab::Copy(
        statePert[lev], ldataNLs_p->nlState, 0, 0, 2, m_nGrowState);
      MultiFab::Saxpy(statePert[lev], delta_pert, *a_v[lev], 0, 0, 2, 0);
    }

    int update_scaling = 0;
    int update_precond = 0;
    nonLinearResidual(
      dtsub, GetVecOfPtrs(statePert), GetVecOfPtrs(residPert), update_scaling,
      update_precond);

    for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
      MultiFab::LinComb(
        *a_Jv[lev], 1.0, residPert[lev], 0, -1.0, ldataNLs_p->nlResid, 0, 0, 2,
        0);
      a_Jv[lev]->mult(-1.0 / delta_pert);
    }
  } else if (m_ef_diffT_jfnk == 2) {
    Vector<MultiFab> statePertPls(finest_level + 1);
    Vector<MultiFab> residPertPls(finest_level + 1);
    Vector<MultiFab> statePertMns(finest_level + 1);
    Vector<MultiFab> residPertMns(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
      statePertPls[lev].define(
        grids[lev], dmap[lev], 2, m_nGrowState, MFInfo(), Factory(lev));
      residPertPls[lev].define(
        grids[lev], dmap[lev], 2, m_nGrowState, MFInfo(), Factory(lev));
      statePertMns[lev].define(
        grids[lev], dmap[lev], 2, m_nGrowState, MFInfo(), Factory(lev));
      residPertMns[lev].define(
        grids[lev], dmap[lev], 2, m_nGrowState, MFInfo(), Factory(lev));
      MultiFab::Copy(
        statePertPls[lev], ldataNLs_p->nlState, 0, 0, 2, m_nGrowState);
      MultiFab::Copy(
        statePertMns[lev], ldataNLs_p->nlState, 0, 0, 2, m_nGrowState);
      MultiFab::Saxpy(statePertPls[lev], delta_pert, *a_v[lev], 0, 0, 2, 0);
      MultiFab::Saxpy(statePertMns[lev], -delta_pert, *a_v[lev], 0, 0, 2, 0);
    }

    int update_scaling = 0;
    int update_precond = 0;
    nonLinearResidual(
      dtsub, GetVecOfPtrs(statePertPls), GetVecOfPtrs(residPertPls),
      update_scaling, update_precond);
    nonLinearResidual(
      dtsub, GetVecOfPtrs(statePertMns), GetVecOfPtrs(residPertMns),
      update_scaling, update_precond);

    for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
      MultiFab::LinComb(
        *a_Jv[lev], 1.0, residPertPls[lev], 0, -1.0, residPertMns[lev], 0, 0, 2,
        0);
      a_Jv[lev]->mult(-0.5 / delta_pert);
    }
  }
}

void
PeleLM::applyPrecond(
  const Vector<MultiFab*>& a_v, const Vector<MultiFab*>& a_Pv)
{
  BL_PROFILE("PeleLMeX::applyPrecond()");

  // Setup aliases and temps
  Vector<MultiFab> nE_al;
  Vector<MultiFab> phiV_al;
  Vector<MultiFab> PnE_al;
  Vector<MultiFab> PphiV_al;
  Vector<std::unique_ptr<MultiFab>> temp(finest_level + 1);
  Vector<std::unique_ptr<MultiFab>> temp2(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    // Aliases
    nE_al.emplace_back(*a_v[lev], amrex::make_alias, 0, 1);
    phiV_al.emplace_back(*a_v[lev], amrex::make_alias, 1, 1);
    PnE_al.emplace_back(*a_Pv[lev], amrex::make_alias, 0, 1);
    PphiV_al.emplace_back(*a_Pv[lev], amrex::make_alias, 1, 1);
    PnE_al[lev].setVal(0.0, 0, 1, 1);
    PphiV_al[lev].setVal(0.0, 0, 1, 1);

    // Temporary data holder
    temp[lev].reset(new MultiFab(
      grids[lev], dmap[lev], 1, nE_al[lev].nGrow(), MFInfo(), Factory(lev)));
    temp[lev]->setVal(0.0);
    temp2[lev].reset(new MultiFab(
      grids[lev], dmap[lev], 1, nE_al[lev].nGrow(), MFInfo(), Factory(lev)));
    temp2[lev]->setVal(0.0);
  }

  // Set LinOps Level BCs
  for (int lev = 0; lev <= finest_level; ++lev) {
    getPrecondOp()->setDiffOpLevelBC(lev, &PnE_al[lev]);
    getPrecondOp()->setDriftOpLevelBC(lev, &PphiV_al[lev]);
    getPrecondOp()->setStildaOpLevelBC(lev, &PphiV_al[lev]);
  }

  // Most inner matrix
  Real S_tol = m_ef_PC_MG_Tol;
  Real S_tol_abs = MLNorm0(GetVecOfConstPtrs(nE_al)) * m_ef_PC_MG_Tol;
  getPrecondOp()->diffOpSolve(
    GetVecOfPtrs(PnE_al), GetVecOfConstPtrs(nE_al), S_tol, S_tol_abs);
  for (int lev = 0; lev <= finest_level; ++lev) {
    MultiFab::Copy(PphiV_al[lev], phiV_al[lev], 0, 0, 1, 0);
  }

  // Pivot second matrix
  for (int lev = 0; lev <= finest_level; ++lev) {
    MultiFab::Saxpy(
      PphiV_al[lev], nE_scale / FphiV_scale, PnE_al[lev], 0, 0, 1, 0);
  }

  // PhiV estimate matrix
  // First scale PphiV
  for (int lev = 0; lev <= finest_level; ++lev) {
    PphiV_al[lev].mult(FphiV_scale / phiV_scale);
  }
  S_tol_abs = MLNorm0(GetVecOfConstPtrs(PphiV_al)) * m_ef_PC_MG_Tol;
  getPrecondOp()->StildaOpSolve(
    GetVecOfPtrs(temp), GetVecOfConstPtrs(PphiV_al), S_tol, S_tol_abs);
  for (int lev = 0; lev <= finest_level; ++lev) {
    MultiFab::Copy(PphiV_al[lev], *temp[lev], 0, 0, 1, 0);
  }

  // Final matrix
  getPrecondOp()->driftOpApply(GetVecOfPtrs(temp), GetVecOfPtrs(PphiV_al));
  S_tol_abs = MLNorm0(GetVecOfConstPtrs(temp)) * m_ef_PC_MG_Tol;
  getPrecondOp()->diffOpSolve(
    GetVecOfPtrs(temp2), GetVecOfConstPtrs(temp), S_tol, S_tol_abs);
  for (int lev = 0; lev <= finest_level; ++lev) {
    MultiFab::Saxpy(PnE_al[lev], -1.0, *temp2[lev], 0, 0, 1, 0);
  }
}

void
PeleLM::nlSolveNorm(const Vector<MultiFab*>& a_MF, Real& r)
{
  r = 0.0;
  int nComp = a_MF[0]->nComp();
  for (int comp = 0; comp < nComp; comp++) {
    Real norm = 0.0;
    for (int lev = 0; lev < a_MF.size(); ++lev) {
      // TODO : norm not weighted by cell size, should it ?
      if (lev != a_MF.size() - 1) {
        norm += MultiFab::Dot(
          *m_coveredMask[lev], *a_MF[lev], comp, *a_MF[lev], comp, 1, 0);
      } else {
        norm += MultiFab::Dot(*a_MF[lev], comp, *a_MF[lev], comp, 1, 0);
      }
    }
    r += norm;
  }
  r = std::sqrt(r);
}
