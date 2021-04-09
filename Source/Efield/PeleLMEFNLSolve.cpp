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

   // Get cell-centered charged species transport coefficients
   calcEFTransport(AmrNewTime);

   // Substepping of non-linear solve
   dtsub = a_dt/ef_substep;

   // Pass t^{n} nE/PhiV from leveldata to leveldatanlsolve
   // t^{n} have been fillpatched already
   for (int lev = 0; lev <= finest_level; ++lev) {
      // Get t^{n} data pointer
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
      // Get nl solve data pointer
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);
      int nGrowNL = 1;
      MultiFab::Copy(ldataNLs_p->nlState, ldata_p->nE, 0, 0, 1, nGrowNL);
      MultiFab::Copy(ldataNLs_p->nlState, ldata_p->phiV, 0, 1, 1, nGrowNL);
   }

   // Gradient of PhiV at t^{n}
   int do_avgDown = 0;     // TODO or should I ?
   auto bcRecPhiV = fetchBCRecArray(PHIV,1);
   getDiffusionOp()->computeGradient(getNLgradPhiVVect(),
                                     {},           // don't need the laplacian out
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
      // TODO scale nl state vector
      if ( ef_verbose ) {
         amrex::Print() << "(" << sstep << ") ne scaling: " << nE_scale << "\n";
         amrex::Print() << "(" << sstep << ") phiV scaling: " << phiV_scale << "\n";
      }    

      // Compute the background charge distribution at curtime
      computeBGcharge(curtime, diffData, advData);

      // Newton initial guess
      // TODO

      // Initial NL residual: update residual scaling and preconditioner
      int update_scaling = 1;
      int update_precond = 1;
      nonLinearResidual(dtsub, getNLstateVect(), getNLresidVect(), update_scaling, update_precond);
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

         // Solve for Newton direction

         // Linesearch & update state

         // Exit condition
         exit_newton = (NK_ite > m_ef_maxNewtonIter);

      } while( !exit_newton );
      NK_tot_count += NK_ite;

      // -----------------
      // Post-Newton
      // Increment the forcing term
      // Unscale nl_state and if not last subcycle update 'old' state
   }

   // Update the state
   //
   Abort();
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
                               const Vector<MultiFab*> &a_nlstate,
                               const Vector<MultiFab*> &a_nlresid,
                               int updateScaling,
                               int updatePrecond)
{
   // Get unscaled copy of the NL state
   int nGrowNLstate = 1;
   Vector<MultiFab> nE(finest_level+1);
   Vector<MultiFab> phiV(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      nE[lev].define(grids[lev],dmap[lev],1,nGrowNLstate,MFInfo(),Factory(lev));
      phiV[lev].define(grids[lev],dmap[lev],1,nGrowNLstate,MFInfo(),Factory(lev));
      MultiFab::Copy(nE[lev],*a_nlstate[lev],0,0,1,nGrowNLstate);
      // TODO nE[lev].mult(nE_scale,0,1,nGrowNLstate);
      MultiFab::Copy(phiV[lev],*a_nlstate[lev],1,0,1,nGrowNLstate);
      // TODO phiV[lev].mult(phiV_scale,0,1,nGrowNLstate);
   }

   // Get L(phiV) and Grad(phiV)
   Vector<MultiFab> laplacian(finest_level+1);
   Vector<Array<MultiFab,AMREX_SPACEDIM>> gradPhiVCur(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      laplacian[lev].define(grids[lev],dmap[lev],1,0,MFInfo(),Factory(lev));
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         const auto& fba = amrex::convert(grids[lev],IntVect::TheDimensionVector(idim));
         gradPhiVCur[lev][idim].define(fba,dmap[lev],1,0,MFInfo(),Factory(lev));
      }
   }
   int do_avgDown = 0;  // TODO or shoud I ?
   auto bcRecPhiV = fetchBCRecArray(PHIV,1);
   getDiffusionOp()->computeGradient(GetVecOfArrOfPtrs(gradPhiVCur),
                                     GetVecOfPtrs(laplacian),
                                     GetVecOfConstPtrs(phiV),
                                     bcRecPhiV[0], do_avgDown);

   VisMF::Write(laplacian[0],"lapPhiV");
   VisMF::Write(gradPhiVCur[0][1],"gradPhiVY");

   // Get nE diffusion term
   Vector<MultiFab> diffnE(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      diffnE[lev].define(grids[lev],dmap[lev],1,0,MFInfo(),Factory(lev));
   }
   auto bcRecnE = fetchBCRecArray(NE,1);
   getDiffusionOp()->computeDiffLap(GetVecOfPtrs(diffnE), 0,
                                    GetVecOfConstPtrs(nE), 0,
                                    {},
                                    GetVecOfConstPtrs(getnEDiffusivityVect(AmrNewTime)), 0,
                                    bcRecnE, 1);

   WriteDebugPlotFile(GetVecOfConstPtrs(diffnE),"diffnE");
   WriteDebugPlotFile(GetVecOfConstPtrs(nE),"nEforadv");

   // Get nE advection term
   // TODO: at some point, switch to using Godunov/MOL
   Vector<MultiFab> advnE(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      advnE[lev].define(grids[lev],dmap[lev],1,0,MFInfo(),Factory(lev));
   }
   getAdvectionTerm(GetVecOfConstPtrs(nE),
                    GetVecOfPtrs(advnE),
                    GetVecOfArrOfConstPtrs(gradPhiVCur));
   WriteDebugPlotFile(GetVecOfConstPtrs(advnE),"advnE");

   // Assemble non-linear residual
   // res(ne(:)) = dt * ( diff(:) + conv(:) + I_R(:) ) - ( ne(:) - ne_old(:) )
   // res(phiv(:)) = \Sum z_k * \tilde Y_k / q_e - ne + Lapl_PhiV
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);  // New time
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);  // Old time
      auto ldataR_p = getLevelDataReactPtr(lev);       // Reaction

      // Get nl solve data pointer
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNLs_p->nlResid,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& I_R_nE   = ldataR_p->I_RnE.const_array(mfi);
         auto const& lapPhiV  = laplacian[lev].const_array(mfi);
         auto const& ne_diff  = diffnE[lev].const_array(mfi);
         auto const& ne_adv   = advnE[lev].const_array(mfi);
         auto const& ne_curr  = nE[lev].const_array(mfi);
         auto const& ne_old   = ldataOld_p->nE.const_array(mfi);
         auto const& charge   = ldataNLs_p->backgroundCharge.const_array(mfi);
         auto const& res_nE   = a_nlresid[lev]->array(mfi,0);
         auto const& res_phiV = a_nlresid[lev]->array(mfi,1);
         Real scalLap         = eps0 * epsr / elemCharge;
         amrex::ParallelFor(bx, [ne_curr,ne_old,lapPhiV,I_R_nE,ne_diff,ne_adv,charge,res_nE,res_phiV,a_dt,scalLap]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            res_nE(i,j,k) = ne_old(i,j,k) - ne_curr(i,j,k) + a_dt * ( ne_diff(i,j,k) + ne_adv(i,j,k) /* TODO + I_R_nE(i,j,k)*/ );
            res_phiV(i,j,k) = lapPhiV(i,j,k) * scalLap - ne_curr(i,j,k) + charge(i,j,k);
         });
      }
   }
   WriteDebugPlotFile(GetVecOfConstPtrs(a_nlresid),"nlResid");

   Abort();
}

void PeleLM::getAdvectionTerm(const Vector<const MultiFab*> &a_nE,
                              const Vector<MultiFab*> &a_advTerm,
                              const Vector<Array<const MultiFab*,AMREX_SPACEDIM>> &a_gPhiVCur)
{

   // Get advection fluxes on all levels
   int nGrow = 0;
   Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      const auto& ba = grids[lev];
      const auto& factory = Factory(lev);
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         fluxes[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                  dmap[lev], 1, nGrow, MFInfo(), factory);
      }
   }

   // nE BCRec
   auto bcRecnE = fetchBCRecArray(NE,1);

   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get t^{n} data pointer
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

      // Get nl solve data pointer
      auto ldataNLs_p = getLevelDataNLSolvePtr(lev);

      // Get the face centered electron mobility
      Array<MultiFab,AMREX_SPACEDIM> mobE_ec = getDiffusivity(lev, 0, 1, bcRecnE, ldata_p->mobE_cc);

      // Get the electron effective velocity
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldataNLs_p->uEffnE[idim],TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            auto const& ueff    = ldataNLs_p->uEffnE[idim].array(mfi);
            auto const& gphi_c  = a_gPhiVCur[lev][idim]->const_array(mfi);
            auto const& gphi_o  = ldataNLs_p->gPhiVOld[idim].const_array(mfi);
            auto const& kappa_e = mobE_ec[idim].const_array(mfi);
            amrex::ParallelFor(bx, [ueff, gphi_c, gphi_o, kappa_e]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               ueff(i,j,k) = /* TODO umac */ - kappa_e(i,j,k) * 0.5 * ( gphi_c(i,j,k) + gphi_o(i,j,k) );
            });
         }
      }
   }
   //WriteDebugPlotFile(GetVecOfConstPtrs(getUeffVect()),"Ueff");

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
      getAdvectionFluxes(lev, GetArrOfPtrs(fluxes[lev]), *a_nE[lev], GetArrOfConstPtrs(ldataNLs_p->uEffnE), bcRecnE[0]);
   }

   // Average down the fluxes
   /*
	for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
      EB_average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                            GetArrOfPtrs(fluxes[lev-1]),
                            refRatio(lev-1),geom[lev-1]);
#else
      average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                         GetArrOfPtrs(fluxes[lev-1]),
                         refRatio(lev-1),geom[lev-1]);
#endif
	}
   */

   // Compute divergence
   int intensiveFluxes = 0;
	fluxDivergence(a_advTerm,0,GetVecOfArrOfPtrs(fluxes),0,1,intensiveFluxes,1.0);
}

void PeleLM::getAdvectionFluxes(int lev,
                                const Array<MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                                const MultiFab &a_nE,
                                const Array<const MultiFab*,AMREX_SPACEDIM> &a_ueff,
                                BCRec bcrec)
{
   const Box& domain = geom[lev].Domain();

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      FArrayBox edgstate[AMREX_SPACEDIM];
      for (MFIter mfi(a_nE,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         const Box& bx  = mfi.tilebox();
         const Box& gbx = mfi.growntilebox(1);
         AMREX_D_TERM( const Box& xbx = mfi.grownnodaltilebox(0,0);,
                       const Box& ybx = mfi.grownnodaltilebox(1,0);,
                       const Box& zbx = mfi.grownnodaltilebox(2,0));

         // data arrays
         auto const& ne_arr = a_nE.const_array(mfi);
         AMREX_D_TERM( Array4<Real> xflux = a_fluxes[0]->array(mfi);,
                       Array4<Real> yflux = a_fluxes[1]->array(mfi);,
                       Array4<Real> zflux = a_fluxes[2]->array(mfi));
         AMREX_D_TERM( Array4<Real const> u = a_ueff[0]->const_array(mfi);,
                       Array4<Real const> v = a_ueff[1]->const_array(mfi);,
                       Array4<Real const> w = a_ueff[2]->const_array(mfi););
         AMREX_D_TERM( edgstate[0].resize(xbx,1);,
                       edgstate[1].resize(ybx,1);,
                       edgstate[2].resize(zbx,1));
         AMREX_D_TERM( Array4<Real> xstate = edgstate[0].array();,
                       Array4<Real> ystate = edgstate[1].array();,
                       Array4<Real> zstate = edgstate[2].array());

         // Predict edge states
         // X
         {
            // BCs
            const Box& edomain = surroundingNodes(domain,0);
            const auto bc_lo = bcrec.lo(0);
            const auto bc_hi = bcrec.hi(0);

            amrex::ParallelFor(xbx, [ne_arr,u,xstate,bc_lo,bc_hi,edomain]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               int idx[3] = {i,j,k};
               bool on_lo = ( ( bc_lo == BCType::ext_dir ) && ( idx[0] <= edomain.smallEnd(0) ) );
               bool on_hi = ( ( bc_hi == BCType::ext_dir ) && ( idx[0] >= edomain.bigEnd(0) ) );
               xstate(i,j,k) = ef_edge_state_extdir(i,j,k,0,on_lo,on_hi,ne_arr,u);
            });
         }
#if (AMREX_SPACEDIM > 1)
         // Y
         {
            // BCs
            const Box& edomain = surroundingNodes(domain,1);
            const auto bc_lo = bcrec.lo(1);
            const auto bc_hi = bcrec.hi(1);

            amrex::ParallelFor(ybx, [ne_arr,v,ystate,bc_lo,bc_hi,edomain]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               int idx[3] = {i,j,k};
               bool on_lo = ( ( bc_lo == BCType::ext_dir ) && ( idx[1] <= edomain.smallEnd(1) ) );
               bool on_hi = ( ( bc_hi == BCType::ext_dir ) && ( idx[1] >= edomain.bigEnd(1) ) );
               ystate(i,j,k) = ef_edge_state_extdir(i,j,k,1,on_lo,on_hi,ne_arr,v);
            });
         }
#if ( AMREX_SPACEDIM ==3 )
         // Z
         {
            // BCs
            const Box& edomain = surroundingNodes(domain,2);
            const auto bc_lo = bcrec.lo(2);
            const auto bc_hi = bcrec.hi(2);

            amrex::ParallelFor(zbx, [ne_arr,w,zstate,bc_lo,bc_hi,edomain]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               int idx[3] = {i,j,k};
               bool on_lo = ( ( bc_lo == BCType::ext_dir ) && ( idx[2] <= edomain.smallEnd(2) ) );
               bool on_hi = ( ( bc_hi == BCType::ext_dir ) && ( idx[2] >= edomain.bigEnd(2) ) );
               zstate(i,j,k) = ef_edge_state_extdir(i,j,k,2,on_lo,on_hi,ne_arr,w);
            });
         }
#endif
#endif

         // Computing fluxes
         amrex::ParallelFor(xbx, [u,xstate,xflux]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {    
            xflux(i,j,k) = u(i,j,k) * xstate(i,j,k);
         });  
#if (AMREX_SPACEDIM > 1)
         amrex::ParallelFor(ybx, [v,ystate,yflux]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {    
            yflux(i,j,k) = v(i,j,k) * ystate(i,j,k);
         });  
#if ( AMREX_SPACEDIM ==3 )
         amrex::ParallelFor(zbx, [w,zstate,zflux]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {    
            zflux(i,j,k) = w(i,j,k) * zstate(i,j,k);
         });  
#endif
#endif
      }
   }
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
   Vector<MultiFab> statePert(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      statePert[lev].define(grids[lev],dmap[lev],2,0,MFInfo(),Factory(lev));
   }

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
      for (int lev = 0; lev < a_MF.size(); ++lev) {
         // TODO : norm not weighted by cell size, should it ?
         norm += MultiFab::Dot(*a_MF[lev],comp,*a_MF[lev],comp,1,0);
      }
      r += norm;
   }
   r = std::sqrt(r);
}
