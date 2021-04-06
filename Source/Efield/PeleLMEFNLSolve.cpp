#include <PeleLM.H>
#include <PeleLM_K.H>
#include <PeleLMEF_K.H>
#include <MLGMRES.H>
#include <DiffusionOp.H>

using namespace amrex;

void PeleLM::implicitNonLinearSolve(int sdcIter,
                                    const Real &a_dt,
                                    std::unique_ptr<AdvanceDiffData> &diffData,
                                    std::unique_ptr<AdvanceAdvData> &advData)
{
   BL_PROFILE_VAR("PeleLM::implicitNonLinearSolve()", implicitNonLinearSolve);

   const Real strt_time = ParallelDescriptor::second();

   //------------------------------------------------------------------------
   // Pre-solve

   // Substepping of non-linear solve
   dtsub = a_dt/ef_substep;

   // Pass t^{n} nE/PhiV from leveldata to leveldatanlsolve
   // t^{n} have been fillpatched already
   // TODO will use nE and phiV at AmrNewTime as container for NL solve state
   //for (int lev = 0; lev <= finest_level; ++lev) {
   //   // Get t^{n} data pointer
   //   auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
   //   // Get t^{n} nl solve data pointer
   //   auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
   //   int nGrowNL = 1;
   //   MultiFab::Copy(ldataNLs_p->nlState, ldata_p->nE, 0, 0, 1, nGrowNL);
   //   MultiFab::Copy(ldataNLs_p->nlState, ldata_p->phiV, 0, 1, 1, nGrowNL);
   //}

   // Gradient of PhiV at t^{n}
   int do_avgDown = 0;
   auto bcRecPhiV = fetchBCRecArray(PHIV,1);
   getDiffusionOp()->computeGradient(getNLgradPhiVVect(),
                                     GetVecOfConstPtrs(getPhiVVect(AmrOldTime)),
                                     bcRecPhiV[0], do_avgDown);

   // Setup MLGMRES
   MLGMRESSolver gmres;
   int GMRES_tot_count = 0;
   if ( !m_ef_use_PETSC_direct ) {
      gmres.define(this,m_ef_GMRES_size,2,1);
      MLJtimesVFunc jtv = &PeleLM::jTimesV;
      gmres.setJtimesV(jtv);
      MLNormFunc normF = &PeleLM::nlSolveNorm;
      gmres.setNorm(normF);
      MLPrecondFunc prec = &PeleLM::applyPrecond;
      gmres.setPrecond(prec);
      gmres.setVerbose(m_ef_GMRES_verbose);
      gmres.setMaxRestart(m_ef_GMRES_maxRst);
   }

   //------------------------------------------------------------------------
   // Outer subcycling loop
   int NK_tot_count = 0;
   for (int sstep = 0; sstep < ef_substep; sstep++) {

      curtime = getTime(0,AmrOldTime) + (sstep+1) * dtsub;

      // -----------------
      // Pre-Newton
      // Set up the NL state scaling
      nE_scale = MLNorm0(GetVecOfConstPtrs(getnEVect(AmrNewTime)));
      nE_scale = (nE_scale > 1.0e-12) ? nE_scale : 1.0;
      phiV_scale = MLNorm0(GetVecOfConstPtrs(getPhiVVect(AmrNewTime)));
      phiV_scale = (phiV_scale > 1.0e-12) ? phiV_scale : 1.0;

      // Compute the background charge distribution at curtime
      computeBGcharge(curtime, diffData, advData);

      // Newton initial guess
      // TODO

      // Initial NL residual: update residual scaling and preconditioner
      int update_scaling = 1;
      int update_precond = 1;
      nonLinearResidual(dtsub, update_scaling, update_precond);
      nlSolveNorm(getNLresidVect(),nl_residNorm);

      // Check for convergence
      Real max_nlres = 0.0;

      // -----------------
      // Newton iteration
      int exit_newton = 0;
      int NK_ite = 0;
      do {
         NK_ite += 1;

         // Verbose
         if ( ef_verbose ) {
            amrex::Print() << " Newton it: " << NK_ite << " L2**2 residual: " << 0.5*nl_residNorm*nl_residNorm
                                                       << ". Linf residual: " << max_nlres << "\n";
         }

      } while( !exit_newton );
      NK_tot_count += NK_ite;

      // -----------------
      // Post-Newton
      // Increment the forcing term
      // Unscale nl_state and if not last subcycle update 'old' state
   }

   // Update the state
}

void PeleLM::computeBGcharge(const Real &a_time,
                             std::unique_ptr<AdvanceDiffData> &diffData,
                             std::unique_ptr<AdvanceAdvData> &advData)
{
   // Get integration dt
   Real dt_int = a_time - getTime(0,AmrOldTime);

   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get data pointers
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);  // Old time species
      auto ldataR_p = getLevelDataReactPtr(lev);       // Reaction
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);   // NL data

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNLs_p->backgroundCharge,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rhoYold  = ldata_p->species.const_array(mfi);
         auto const& adv_arr  = advData->AofS[lev].const_array(mfi,FIRSTSPEC);
         auto const& dn_arr   = diffData->Dn[lev].const_array(mfi);
         auto const& dnp1_arr = diffData->Dnp1[lev].const_array(mfi);
         auto const& dhat_arr = diffData->Dhat[lev].const_array(mfi);
         auto const& rhoYdot  = ldataR_p->I_R.const_array(mfi);
         auto const& charge   = ldataNLs_p->backgroundCharge.array(mfi);
         Real        factor = 1.0 / elemCharge;
         amrex::ParallelFor(bx, [dt_int, rhoYold, adv_arr, dn_arr, dnp1_arr, dhat_arr, rhoYdot, charge, factor, zk=zk]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            charge(i,j,k) = 0.0;
            for (int n = 0; n < NUM_SPECIES; n++) {
               Real rhoYprov = rhoYold(i,j,k,n) + dt_int * ( adv_arr(i,j,k,n) +
                                                             0.5 * ( dn_arr(i,j,k,n) - dnp1_arr(i,j,k,n) ) +
                                                             dhat_arr(i,j,k,n) +
                                                             rhoYdot(i,j,k,n) );
               rhoYprov = amrex::max(rhoYprov,0.0);
               charge(i,j,k) += zk[n] * rhoYprov;
            }
            charge(i,j,k) *= factor;
         });
      }
   }
}

void PeleLM::nonLinearResidual(const Real &a_dt,
                              int updateScaling,
                              int updatePrecond)
{
}

void PeleLM::jTimesV(const Vector<MultiFab*> &a_v,
                     const Vector<MultiFab*> &a_Jv)
{
   Real vNorm;
   nlSolveNorm(a_v,vNorm);

   // v is zero, Jv is zero and return
   if ( vNorm == 0.0 ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         a_Jv[lev]->setVal(0.0);
      }
   }

   Real delta_pert = m_ef_lambda_jfnk * ( m_ef_lambda_jfnk + nl_stateNorm / vNorm );

}

void PeleLM::applyPrecond(const Vector<MultiFab*> &a_v,
                          const Vector<MultiFab*> &a_Pv)
{
}

void PeleLM::nlSolveNorm(const Vector<MultiFab*> &a_MF, Real &r)
{
   r = 0.0;
   int nComp = a_MF[0]->nComp();
   for ( int comp = 0; comp < nComp; comp++ ) {
      Real norm = 0.0;
      for (int lev = 0; lev <= a_MF.size(); ++lev) {
         // TODO : norm not weighted by cell size, should it ?
         norm += MultiFab::Dot(*a_MF[lev],comp,*a_MF[lev],comp,1,0);
      }
      r += norm;
   }
   r = std::sqrt(r);
}
