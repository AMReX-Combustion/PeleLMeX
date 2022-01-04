#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

void PeleLM::setThermoPress(const TimeStamp &a_time) {
   BL_PROFILE("PeleLM::setThermoPress()");

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   for (int lev = 0; lev <= finest_level; ++lev) {
      setThermoPress(lev, a_time);
   }

   averageDownRhoRT(a_time);
}

void PeleLM::setThermoPress(int lev, const TimeStamp &a_time) {

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = ldata_p->state.const_array(mfi,DENSITY);
         auto const& rhoY    = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& T       = ldata_p->state.const_array(mfi,TEMP);
         auto const& P       = ldata_p->state.array(mfi,RHORT);

         amrex::ParallelFor(bx, [rho, rhoY, T, P]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getPGivenRTY( i, j, k, rho, rhoY, T, P );
         });
      }
   }
}

void PeleLM::calcDivU(int is_init,
                      int computeDiff,
                      int do_avgDown,
                      const TimeStamp &a_time,
                      std::unique_ptr<AdvanceDiffData> &diffData)
{
   BL_PROFILE("PeleLM::calcDivU()");

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   // If requested, compute diffusion terms
   // otherwise assumes it has already been computed and stored in the proper container of diffData
   if (computeDiff) {
      calcDiffusivity(a_time);
      computeDifferentialDiffusionTerms(a_time, diffData, is_init);
   }

   // Assemble divU on each level
   for (int lev = 0; lev <= finest_level; lev++ ) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

      // Explicitly get old and new data for time derivatives
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

      MultiFab RhoYdot;
      if ( m_do_react && !m_skipInstantRR ) {
         if (is_init) {          // Either pre-divU, divU or press initial iterations
            if ( m_dt > 0.0 ) {  // divU ite   -> use I_R
               auto ldataR_p = getLevelDataReactPtr(lev);
               RhoYdot.define(grids[lev],dmap[lev],nCompIR(),0);
               MultiFab::Copy(RhoYdot,ldataR_p->I_R,0,0,nCompIR(),0);
            } else {             // press ite  -> set to zero
               RhoYdot.define(grids[lev],dmap[lev],nCompIR(),0);
               RhoYdot.setVal(0.0);
            }
         } else {                // Regular    -> use instantaneous RR
            RhoYdot.define(grids[lev],dmap[lev],nCompIR(),0);
#ifdef PELE_USE_EFIELD
            computeInstantaneousReactionRateEF(lev, a_time, &RhoYdot);
#else
            computeInstantaneousReactionRate(lev, a_time, &RhoYdot);
#endif
         }
      }

      //----------------------------------------------------------------
#ifdef AMREX_USE_EB
      // Get EBFact
      const auto& ebfact = EBFactory(lev);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->divu, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
         auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
         auto const& flag    = flagfab.const_array();
#endif

         auto const& rhoY     = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& T        = ldata_p->state.const_array(mfi,TEMP);
         auto const& SpecD    = ( a_time == AmrOldTime ) ? diffData->Dn[lev].const_array(mfi,0)
                                                         : diffData->Dnp1[lev].const_array(mfi,0);
         auto const& Fourier  = ( a_time == AmrOldTime ) ? diffData->Dn[lev].const_array(mfi,NUM_SPECIES)
                                                         : diffData->Dnp1[lev].const_array(mfi,NUM_SPECIES);
         auto const& DiffDiff = ( a_time == AmrOldTime ) ? diffData->Dn[lev].const_array(mfi,NUM_SPECIES+1)
                                                         : diffData->Dnp1[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& r        = (m_do_react && !m_skipInstantRR) ? RhoYdot.const_array(mfi)
                                             : ldata_p->state.const_array(mfi,FIRSTSPEC);   // Dummy unused Array4
         auto const& divu     = ldata_p->divu.array(mfi);
         int use_react        = (m_do_react && !m_skipInstantRR) ? 1 : 0;
         amrex::Real massFluxBalance = AMREX_D_TERM(  m_domainMassFlux[0] + m_domainMassFlux[1],
                                       + m_domainMassFlux[2] + m_domainMassFlux[3],               
                                       + m_domainMassFlux[4] + m_domainMassFlux[5]);
         amrex::Real rhoHFluxBalance = AMREX_D_TERM(  m_domainRhoHFlux[0] + m_domainRhoHFlux[1],
                                       + m_domainRhoHFlux[2] + m_domainRhoHFlux[3],               
                                       + m_domainRhoHFlux[4] + m_domainRhoHFlux[5]);
         Real rhoYFluxBalance[NUM_SPECIES] = {0.0};
         for (int n = 0; n < NUM_SPECIES; n++){
            rhoYFluxBalance[n] = AMREX_D_TERM(  m_domainRhoYFlux[n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[1+n*2*AMREX_SPACEDIM],
                                              + m_domainRhoYFlux[2+n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[3+n*2*AMREX_SPACEDIM],
                                              + m_domainRhoYFlux[4+n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[5+n*2*AMREX_SPACEDIM]); 
         }

         amrex::Real massFluxBalance = AMREX_D_TERM(  m_domainMassFlux[0] + m_domainMassFlux[1],
                                       + m_domainMassFlux[2] + m_domainMassFlux[3],               
                                       + m_domainMassFlux[4] + m_domainMassFlux[5]);
         amrex::Real rhoHFluxBalance = AMREX_D_TERM(  m_domainRhoHFlux[0] + m_domainRhoHFlux[1],
                                       + m_domainRhoHFlux[2] + m_domainRhoHFlux[3],               
                                       + m_domainRhoHFlux[4] + m_domainRhoHFlux[5]);
         Real rhoYFluxBalance[NUM_SPECIES] = {0.0};
         for (int n = 0; n < NUM_SPECIES; n++){
            rhoYFluxBalance[n] = AMREX_D_TERM(  m_domainRhoYFlux[n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[1+n*2*AMREX_SPACEDIM],
                                              + m_domainRhoYFlux[2+n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[3+n*2*AMREX_SPACEDIM],
                                              + m_domainRhoYFlux[4+n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[5+n*2*AMREX_SPACEDIM]);            
         }

#ifdef AMREX_USE_EB
         if (flagfab.getType(bx) == FabType::covered) {             // Covered boxes
             amrex::ParallelFor(bx, [divu]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 divu(i,j,k) = 0.0;
             });
         } else if (flagfab.getType(bx) != FabType::regular ) {     // EB containing boxes 
             amrex::ParallelFor(bx, [ rhoY, T, SpecD, Fourier, DiffDiff, r, rhoHFluxBalance, rhoYFluxBalance, divu, use_react, flag]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                if ( flag(i,j,k).isCovered() ) {
                    divu(i,j,k) = 0.0;
                } else {
                    compute_divu( i, j, k, rhoY, T, SpecD, Fourier, DiffDiff, r, rhoHFluxBalance, rhoYFluxBalance, divu, use_react );
                }
             });
         } else
#endif
         {
             amrex::ParallelFor(bx, [ rhoY, T, SpecD, Fourier, DiffDiff, r, rhoHFluxBalance, rhoYFluxBalance, divu, use_react]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                compute_divu( i, j, k, rhoY, T, SpecD, Fourier, DiffDiff, r, rhoHFluxBalance, rhoYFluxBalance, divu, use_react );
             });
         }
      }
   }

   // Average down divU
   if ( do_avgDown ) {
      for (int lev = finest_level; lev > 0; --lev) {
         auto ldataFine_p = getLevelDataPtr(lev,a_time);
         auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
         EB_average_down(ldataFine_p->divu,
                         ldataCrse_p->divu,
                         0,1,refRatio(lev-1));
#else
         average_down(ldataFine_p->divu,
                      ldataCrse_p->divu,
                      0,1,refRatio(lev-1));
#endif
      }
   }

   // fillPatch a_time divu to get properly filled ghost cells
   int nGrowDivu = 1;   // TODO: need to make sure it's consistent across
   for (int lev = 0; lev <= finest_level; ++lev) {
      Real time = getTime(lev,a_time);
      auto ldata_p = getLevelDataPtr(lev,a_time);
      fillpatch_divu(lev,time,ldata_p->divu,nGrowDivu);
   }
}

void PeleLM::setTemperature(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::setTemperature()", setTemperature);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   for (int lev = 0; lev <= finest_level; ++lev) {
      setTemperature(lev, a_time);
   }
}

void PeleLM::setTemperature(int lev, const TimeStamp &a_time) {

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = ldata_p->state.const_array(mfi,DENSITY);
         auto const& rhoY    = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& rhoh    = ldata_p->state.const_array(mfi,RHOH);
         auto const& T       = ldata_p->state.array(mfi,TEMP);

         amrex::ParallelFor(bx, [rho, rhoY, rhoh, T]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getTfromHY( i, j, k, rho, rhoY, rhoh, T);
         });
      }
   }
}

void PeleLM::calc_dPdt(const TimeStamp &a_time,
                       const Vector<MultiFab*> &a_dPdt)
{
   BL_PROFILE("PeleLM::calc_dPdt()");

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   for (int lev = 0; lev <= finest_level; ++lev) {
      calc_dPdt(lev, a_time, a_dPdt[lev]);
#ifdef AMRE_USE_EB
      EB_set_covered(a_dPdt[lev],0.0);
#endif
   }

   // Fill ghost cell(s)
   if (a_dPdt[0]->nGrow() > 0) {
      fillpatch_forces(m_cur_time,a_dPdt,a_dPdt[0]->nGrow());
   }
}

void PeleLM::calc_dPdt(int lev,
                       const TimeStamp &a_time,
                       MultiFab* a_dPdt)
{

   auto ldata_p = getLevelDataPtr(lev,a_time);

   // Use new ambient pressure to compute dPdt
   Real p_amb = m_pNew;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(*a_dPdt,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& dPdt  = a_dPdt->array(mfi);
      auto const& P     = ldata_p->state.const_array(mfi,RHORT);
      amrex::ParallelFor(bx, [dPdt, P, p_amb, dt=m_dt, dpdt_fac=m_dpdtFactor]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         dPdt(i,j,k) = (P(i,j,k) - p_amb) / ( dt * P(i,j,k) ) * dpdt_fac;
      });
   }
}

Real
PeleLM::adjustPandDivU(std::unique_ptr<AdvanceAdvData> &advData)
{
    BL_PROFILE("PeleLM::adjustPandDivU()");

    Vector<std::unique_ptr<MultiFab>> ThetaHalft(finest_level+1);

    // Get theta = 1 / (\Gamma * P_amb) at half time
    for (int lev = 0; lev <= finest_level; ++lev) {

        auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
        auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

        ThetaHalft[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0, MFInfo(), *m_factory[lev]));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(*ThetaHalft[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& rhoYo  = ldataOld_p->state.const_array(mfi,FIRSTSPEC);
            auto const& rhoYn  = ldataNew_p->state.const_array(mfi,FIRSTSPEC);
            auto const& T_o    = ldataOld_p->state.const_array(mfi,TEMP);
            auto const& T_n    = ldataNew_p->state.const_array(mfi,TEMP);
            auto const& theta  = ThetaHalft[lev]->array(mfi);
            amrex::ParallelFor(bx, [=,pOld=m_pOld,pNew=m_pNew]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                Real gammaInv_o = getGammaInv(i,j,k,rhoYo,T_o); 
                Real gammaInv_n = getGammaInv(i,j,k,rhoYn,T_n); 
                theta(i,j,k) = 0.5 * (gammaInv_o/pOld + gammaInv_n/pNew);
            });
        }
    }

      // NOTE: This are probably obsolete and should be removed in the next push...
//     // Add in source terms to divu here...
//     for (int lev = 0; lev <= finest_level; ++lev) {
// #ifdef AMREX_USE_OMP
// #pragma omp parallel if (Gpu::notInLaunchRegion())
// #endif
//         auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
//         auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

//         for (MFIter mfi(*ThetaHalft[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
//         {
//             const Box& bx = mfi.tilebox();
//             auto const& theta    = ThetaHalft[lev]->const_array(mfi);
//             auto const& macdivU  = advData->mac_divu[lev].array(mfi);
//             auto const& rhoY     = ldataNew_p->species.const_array(mfi);
//             auto const& T        = ldataNew_p->temp.const_array(mfi);
//             amrex::Real massFluxBalance = AMREX_D_TERM(  m_domainMassFlux[0] + m_domainMassFlux[1],
//                                        + m_domainMassFlux[2] + m_domainMassFlux[3],               
//                                        + m_domainMassFlux[4] + m_domainMassFlux[5]);
//             // amrex::Print() << "massFluxBalance = " << massFluxBalance << std::endl;
//             amrex::Real rhoHFluxBalance = AMREX_D_TERM(  m_domainRhoHFlux[0] + m_domainRhoHFlux[1],
//                                           + m_domainRhoHFlux[2] + m_domainRhoHFlux[3],               
//                                           + m_domainRhoHFlux[4] + m_domainRhoHFlux[5]);
//             // amrex::Print() << "rhoHFluxBalance = " << rhoHFluxBalance << std::endl;
//             Real rhoYFluxBalance[NUM_SPECIES] = {0.0};
//             for (int n = 0; n < NUM_SPECIES; n++){
//                rhoYFluxBalance[n] = AMREX_D_TERM(  m_domainRhoYFlux[n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[1+n*2*AMREX_SPACEDIM],
//                                                  + m_domainRhoYFlux[2+n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[3+n*2*AMREX_SPACEDIM],
//                                                  + m_domainRhoYFlux[4+n*2*AMREX_SPACEDIM] + m_domainRhoYFlux[5+n*2*AMREX_SPACEDIM]); 
//                // amrex::Print() << "rhoYFluxBalance[" << n << "] = " << rhoYFluxBalance[n] << std::endl;
//             }

//             // amrex::ParallelFor(bx, [=,dp0dt=m_dp0dt, &rhoY, &T, &rhoHFluxBalance, &rhoYFluxBalance]
//             // AMREX_GPU_DEVICE (int i, int j, int k) noexcept
//             // {
//             //    // compute terms for divu source
//             //    auto eos = pele::physics::PhysicsType::eos();
//             //    amrex::Real mwtinv[NUM_SPECIES] = {0.0};
//             //    eos.inv_molecular_weight(mwtinv);
//             //    // Get rho & Y from rhoY
//             //    amrex::Real rho = 0.0_rt;
//             //    for (int n = 0; n < NUM_SPECIES; n++) {
//             //       rho += rhoY(i,j,k,n);
//             //    }
//             //    amrex::Real rhoinv = 1.0_rt / rho;
//             //    amrex::Real y[NUM_SPECIES] = {0.0};
//             //    for (int n = 0; n < NUM_SPECIES; n++) {
//             //       y[n] = rhoY(i,j,k,n) * rhoinv;
//             //    }
//             //    // get cpmix from T and Y
//             //    amrex::Real cpmix = 0.0_rt;
//             //    eos.TY2Cp(T(i,j,k), y, cpmix);
//             //    amrex::Real Wbar = 0.0_rt;
//             //    eos.Y2WBAR(y, Wbar);
//             //    amrex::Real hi[NUM_SPECIES] = {0.0};
//             //    eos.T2Hi(T(i,j,k), hi);
//             //    cpmix *= 1.0e-4_rt;                                // CGS -> MKS conversion
//             //    for (int n = 0; n < NUM_SPECIES; n++) {
//             //       hi[n] *= 1.0e-4_rt;                             // CGS -> MKS conversion
//             //    }
//             //    amrex::Real denominv = 1.0_rt / ( rho * cpmix * T(i,j,k) );
//             //    //Need to use the face rho, T, for the source terms -- different from volumetric sources terms

//             //    // // Add in the source terms
//             //    //  macdivU(i,j,k) += rhoHFluxBalance * denominv;
//             //    //  for (int n = 0; n < NUM_SPECIES; n++) {
//             //    //    macdivU(i,j,k) += rhoYFluxBalance[n] * mwtinv[n] * Wbar * rhoinv;
//             //    //    macdivU(i,j,k) -= (hi[n] * rhoYFluxBalance[n]) * denominv;
//             //    //  }


//             // });
//          }
//       }


    Vector<MultiFab> dummy(finest_level+1);
    for (int lev = 0; lev <= finest_level; ++lev) {
       dummy[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *m_factory[lev]);
       dummy[lev].setVal(1.0);
    }

    m_uncoveredVol = MFSum(GetVecOfConstPtrs(dummy),0);

    // Get the mean mac_divu (Sbar) and mean theta
    Real Sbar = MFSum(GetVecOfConstPtrs(advData->mac_divu),0);
    Sbar /= m_uncoveredVol;
    Real Thetabar = MFSum(GetVecOfConstPtrs(ThetaHalft),0);
    Thetabar /= m_uncoveredVol;

    // Adjust
    for (int lev = 0; lev <= finest_level; ++lev) {
        // ThetaHalft is now delta_theta
        ThetaHalft[lev]->plus(-Thetabar,0,1);
        // mac_divu is now delta_S
        advData->mac_divu[lev].plus(-Sbar,0,1);
    }
  
    // Compute 1/Volume * int(U_inflow)dA
    amrex::Real umacFluxBalance = AMREX_D_TERM(  m_domainUmacFlux[0] + m_domainUmacFlux[1],
                                       + m_domainUmacFlux[2] + m_domainUmacFlux[3],               
                                       + m_domainUmacFlux[4] + m_domainUmacFlux[5]);
    Real divu_vol = umacFluxBalance/m_uncoveredVol;
    // Real divu_vol = 0.0;

    amrex::Print() << "Total mass_old = " << m_massOld << std::endl;
    amrex::Print() << "Pressure old = " << m_pOld << std::endl;
    amrex::Print() << "umacFluxBalance = " << umacFluxBalance << std::endl;
    amrex::Print() << "divu_vol = " << divu_vol << std::endl;

    // Advance the ambient pressure
    m_pNew = m_pOld + m_dt * (Sbar/Thetabar - divu_vol/Thetabar);
    m_dp0dt = Sbar/Thetabar - divu_vol/Thetabar;

    amrex::Print() << "Sbar/Thetabar = " << Sbar/Thetabar << std::endl;
    amrex::Print() << "divu_vol/Thetabar = " << divu_vol/Thetabar << std::endl;

    // subtract \tilde{theta} * Sbar / Thetabar from divu
    for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif

        for (MFIter mfi(*ThetaHalft[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto const& theta    = ThetaHalft[lev]->const_array(mfi);
            auto const& macdivU  = advData->mac_divu[lev].array(mfi);
            amrex::ParallelFor(bx, [=,dp0dt=m_dp0dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               // Do the closed chamber pressure correction
                // macdivU(i,j,k) -= theta(i,j,k) * dp0dt;
               macdivU(i,j,k) -= (theta(i,j,k) * Sbar/Thetabar - divu_vol*(1 + theta(i,j,k)/Thetabar));
            });
        }
    }

    if (m_verbose > 2 ) {
        Print() << " >> Closed chamber pOld: " << m_pOld << ", pNew: " << m_pNew << ", dp0dt: " << m_dp0dt << "\n";
    }

    // Return Sbar so that we'll add it back to mac_divu after the MAC projection
    return Sbar;
}
