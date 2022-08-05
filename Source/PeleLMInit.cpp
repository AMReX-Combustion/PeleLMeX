#include <PeleLM.H>
#include <pelelm_prob.H>
#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#endif

using namespace amrex;

void PeleLM::Init() {
   BL_PROFILE_VAR("PeleLM::Init()", Init);

   // Open temporals file
   openTempFile();

   // Check run parameters
   checkRunParams();

   // Initialize data
   initData();
}

void PeleLM::MakeNewLevelFromScratch( int lev,
                                            amrex::Real time,
                                      const amrex::BoxArray& ba,
                                      const amrex::DistributionMapping& dm) {
   BL_PROFILE_VAR("PeleLM::MakeNewLevelFromScratch()", MakeNewLevelFromScratch);

   if (m_verbose > 0) {
      amrex::Print() << " Making new level " << lev << " from scratch" << std::endl;
      if (m_verbose > 2 && lev > 0) {
         auto const dx = geom[lev].CellSizeArray();
         Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
         amrex::Print() << " with " << ba.numPts() << " cells,"
                        << " over " << ba.numPts() * vol / geom[0].ProbSize() * 100 << "% of the domain \n";
      }
      if (m_verbose > 3 && lev > 0) {
         amrex::Print() << " with BoxArray " << ba << std::endl;
      }
   }

   // Pass Box and Dmap to AmrCore
   SetBoxArray(lev, ba);
   SetDistributionMap(lev, dm);

   // Define the FAB Factory
#ifdef AMREX_USE_EB
   m_factory[lev] = makeEBFabFactory(geom[lev], grids[lev], dmap[lev],
                                     {6,6,6},
                                     EBSupport::full);
#else
   m_factory[lev].reset(new FArrayBoxFactory());
#endif

   // Initialize the LevelData
   m_leveldata_old[lev].reset(new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                            m_incompressible, m_has_divu,
                                            m_nAux, m_nGrowState));
   m_leveldata_new[lev].reset(new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                            m_incompressible, m_has_divu,
                                            m_nAux, m_nGrowState));

   if (max_level > 0 && lev != max_level) {
      m_coveredMask[lev].reset(new iMultiFab(grids[lev], dmap[lev], 1, 0));
      m_resetCoveredMask = 1;
   }
   if (m_do_react) {
      m_leveldatareact[lev].reset(new LevelDataReact(grids[lev], dmap[lev], *m_factory[lev]));
      m_leveldatareact[lev]->functC.setVal(0.0);
   }

#ifdef PELE_USE_EFIELD
   m_leveldatanlsolve[lev].reset(new LevelDataNLSolve(grids[lev], dmap[lev], *m_factory[lev], m_nGrowState));
   if (m_do_extraEFdiags) {
      m_ionsFluxes[lev].reset(new MultiFab(grids[lev], dmap[lev], NUM_IONS*AMREX_SPACEDIM, 0));
   }
#endif

   // Fill the initial solution (if not restarting)
   if (m_restart_chkfile.empty()) {
      if (m_restart_pltfile.empty()) {
          initLevelData(lev);
      } else {
          initLevelDataFromPlt(lev, m_restart_pltfile);
      }
   }

   // Times
   m_t_new[lev] = time;
   m_t_old[lev] = time - 1.0e200;

   // Mac projector
#if AMREX_USE_EB
   macproj.reset(new Hydro::MacProjector(Geom(0,finest_level),
                                         MLMG::Location::FaceCentroid,  // Location of mac velocity
                                         MLMG::Location::FaceCentroid,  // Location of beta
                                         MLMG::Location::CellCenter) ); // Location of solution variable phi
#else
   macproj.reset(new Hydro::MacProjector(Geom(0,finest_level)));
#endif
   m_macProjOldSize = finest_level+1;
   m_extSource[lev].reset(new MultiFab(grids[lev], dmap[lev], NVAR, amrex::max(m_nGrowAdv, m_nGrowMAC),
                                       MFInfo(), *m_factory[lev]));
   m_extSource[lev]->setVal(0.);

#if AMREX_USE_EB
   if ( lev == 0 && m_signDistNeeded) {
      // Set up CC signed distance container to control EB refinement
      m_signedDist0.reset(new MultiFab(grids[lev], dmap[lev], 1, 1, MFInfo(), *m_factory[lev]));

      // Estimate the maximum distance we need in terms of level 0 dx:
      Real extentFactor = static_cast<Real>(nErrorBuf(0));
      for (int ilev = 1; ilev <= max_level; ++ilev) {
          extentFactor += static_cast<Real>(nErrorBuf(ilev)) / std::pow(static_cast<Real>(refRatio(ilev-1)[0]),
                                                                        static_cast<Real>(ilev));
      }
      extentFactor *= std::sqrt(2.0) * m_derefineEBBuffer;  // Account for diagonals

      MultiFab signDist(convert(grids[0],IntVect::TheUnitVector()),dmap[0],1,1,MFInfo(),EBFactory(0));
      FillSignedDistance(signDist,true);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*m_signedDist0,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.growntilebox();
         auto const& sd_cc = m_signedDist0->array(mfi);
         auto const& sd_nd = signDist.const_array(mfi);
         amrex::ParallelFor(bx, [sd_cc, sd_nd]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            amrex::Real fac = AMREX_D_PICK(0.5,0.25,0.125);
            sd_cc(i,j,k) = AMREX_D_TERM(  sd_nd(i,j,k)   + sd_nd(i+1,j,k),
                                        + sd_nd(i,j+1,k) + sd_nd(i+1,j+1,k),
                                        + sd_nd(i,j,k+1) + sd_nd(i+1,j,k+1)
                                        + sd_nd(i,j+1,k+1) + sd_nd(i+1,j+1,k+1));
            sd_cc(i,j,k) *= fac;
         });
      }
      m_signedDist0->FillBoundary(geom[0].periodicity());
      extendSignedDistance(m_signedDist0.get(), extentFactor);
   }
#endif
}

void PeleLM::initData() {
   BL_PROFILE_VAR("PeleLM::initData()", initData);

   if (m_restart_chkfile.empty()) {

      //----------------------------------------------------------------
      if (!m_initial_grid_file.empty()) {
         InitFromGridFile(m_cur_time);
      } else {
         // This is an AmrCore member function which recursively makes new levels
         // with MakeNewLevelFromScratch.
         InitFromScratch(m_cur_time);
      }
      resetCoveredMask();
      updateDiagnostics();

#ifdef PELELM_USE_SPRAY
      if (do_spray_particles) {
        initSprays();
      }
#endif

      //----------------------------------------------------------------
      // Set typical values
      int is_init = 1;
      setTypicalValues(AmrNewTime, is_init);

      // initiliaze temporals
      initTemporals(AmrNewTime);

#ifdef AMREX_USE_EB
      //----------------------------------------------------------------
      // Initial redistribution
      initCoveredState();
      initialRedistribution();
#endif

      //----------------------------------------------------------------
      // AverageDown and FillPatch the NewState
      averageDownState(AmrNewTime);
      fillPatchState(AmrNewTime);

      //----------------------------------------------------------------
      // If performing UnitTest, let's stop here
      if (runMode() != "normal") {
         return;
      }

#ifdef PELE_USE_EFIELD
      poissonSolveEF(AmrNewTime);
      fillPatchPhiV(AmrNewTime);
#endif

      // Post data Init time step estimate
      // TODO : this estimate is probably useless
      Real dtInit = computeDt(is_init,AmrNewTime);
      Print() << " Initial dt: " << dtInit << "\n";
      //WriteDebugPlotFile(GetVecOfConstPtrs(getStateVect(AmrNewTime)),"InitSol");

      // Subcycling IAMR/PeleLM first does a projection with no reaction divU
      // which can make the dt for evaluating I_R better
      if (m_has_divu) {
         int is_initialization = 1;             // Yes we are
         int computeDiffusionTerm = 1;          // Needed here
         int do_avgDown = 1;                    // Always

         // Light version of the diffusion data container
         std::unique_ptr<AdvanceDiffData> diffData;
         diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory,
                        m_nGrowAdv, m_use_wbar, is_initialization));
         calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
      }
      initialProjection();
     
      // If gravity is used, do initial pressure projection to get the hydrostatic pressure
      if (std::abs(m_gravity.sum()) > 0.0) {
         initialPressProjection();
      }

      // Post data Init time step estimate
      m_dt = computeDt(is_init,AmrNewTime);
      Print() << " Initial dt: " << m_dt << "\n";

      //----------------------------------------------------------------
      // Initial velocity projection iterations
      if (m_do_init_proj) {
         for (int iter = 0; iter < m_numDivuIter; iter++) {
            if (m_do_react) {
               // The new level data has been filled above
               // Copy new -> old since old used in advanceChemistry
               copyStateNewToOld();

               for (int lev = finest_level; lev >= 0; --lev) {
                  // Setup empty forcing
                  MultiFab Forcing(grids[lev],dmap[lev],nCompForcing(),0);
                  Forcing.setVal(0.0);

                  if (lev != finest_level) {
                     // On all but the finest level, average down I_R
                     // and use advanceChemistry with chem BoxArray
                     std::unique_ptr<MultiFab> avgDownIR;
                     avgDownIR.reset(new MultiFab(grids[lev],dmap[lev],nCompIR(),0));
                     avgDownIR->setVal(0.0);
                     auto ldataRFine_p   = getLevelDataReactPtr(lev+1);
#ifdef AMREX_USE_EB
                     EB_average_down(ldataRFine_p->I_R,
                                     *avgDownIR,
                                     0,nCompIR(),refRatio(lev));
#else
                     average_down(ldataRFine_p->I_R,
                                  *avgDownIR,
                                  0,nCompIR(),refRatio(lev));
#endif
                     // Call advanceChemistry
                     advanceChemistry(lev, dtInit/2.0, Forcing, avgDownIR.get());
                  } else {
                     // Call advanceChemistry
                     advanceChemistry(lev, dtInit/2.0, Forcing);
                  }
               }
               // Copy back old -> new
               copyStateOldToNew();
            }
            if (m_has_divu) {
               int is_initialization = 1;             // Yes we are
               int computeDiffusionTerm = 1;          // Needed here
               int do_avgDown = 1;                    // Always

               // Light version of the diffusion data container
               std::unique_ptr<AdvanceDiffData> diffData;
               diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory,
                              m_nGrowAdv, m_use_wbar, is_initialization));
               calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
            }

            initialProjection();
         }
         if ( m_numDivuIter == 0 ) {
            for (int lev = 0; lev <= finest_level; ++lev) {
               auto ldataR_p   = getLevelDataReactPtr(lev);
               ldataR_p->I_R.setVal(0.0);
            }
         }
      }

      //----------------------------------------------------------------
      // Initial pressure iterations
      initialIterations();

      m_nstep = 0;

      if (m_do_temporals) {
         writeTemporals();
      }

      if (m_plot_int > 0 || m_plot_per_approx > 0. || m_plot_per_exact > 0.) {
         WritePlotFile();
      }
      if (m_check_int > 0 ) {
         WriteCheckPointFile();
      }

   } else {
      //----------------------------------------------------------------
      // Read starting configuration from chk file.
      ReadCheckPointFile();

#ifdef PELELM_USE_SPRAY
      if (do_spray_particles) {
        sprayRestart();
      }
#endif
#ifdef PELE_USE_EFIELD
      // If restarting from a non efield simulation
      if (m_restart_nonEF) {
         // either pass Y_ne -> nE or initialize nE for electro-neutral
         if (m_restart_electroneutral) {
            initializeElectronNeutral();
         } else {
            initializeElectronFromMassFraction();
         }

         // do an initial Poisson solve
         fillPatchPhiV(AmrNewTime);
         poissonSolveEF(AmrNewTime);

         // Reset time data
         if ( m_restart_resetTime ) {
            m_nstep = 0;
            m_cur_time = 0.0;
            for (int lev = 0; lev <= finest_level; ++lev) {
                m_t_new[lev] = 0.0;
                m_t_old[lev] = -1.0e200;
            }
            m_dt = -1.0;
            int is_init = 1;
            Real dtInit = computeDt(is_init,AmrNewTime);
            Print() << " Initial dt: " << dtInit << "\n";
         }

         // Let's write the initial condition
         if (m_plot_int > 0 ) {
            WritePlotFile();
         }
      }
#endif

      // Generate the covered cell mask
      m_resetCoveredMask = 1;
      resetCoveredMask();
      updateDiagnostics();
   }

}

void PeleLM::initLevelData(int lev) {
   BL_PROFILE_VAR("PeleLM::initLevelData()", initLevelData);

   // Get level data
   auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
   const auto geomdata = geom[lev].data();

   // Pressure and pressure gradients to zero
   ldata_p->press.setVal(0.0);
   ldata_p->gp.setVal(0.0);

   // Prob/PMF datas
   ProbParm const* lprobparm = prob_parm_d;
   pele::physics::PMF::PmfData::DataContainer const* lpmfdata   = pmf_data.getDeviceData();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      FArrayBox DummyFab(bx,1);
      auto  const &state_arr   = ldata_p->state.array(mfi);
      auto  const &aux_arr   = (m_nAux > 0) ? ldata_p->auxiliaries.array(mfi) : DummyFab.array();
      amrex::ParallelFor(bx, [=,m_incompressible=m_incompressible]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         pelelm_initdata(i, j, k, m_incompressible, state_arr, aux_arr,
                         geomdata, *lprobparm, lpmfdata);
      });
   }

   if (!m_incompressible) {
      // Initialize thermodynamic pressure
      setThermoPress(lev, AmrNewTime);
      if (m_has_divu) {
         ldata_p->divu.setVal(0.0);
      }
   }
}

void PeleLM::initialIterations() {
   BL_PROFILE_VAR("PeleLM::initialIterations()", initialIterations);

   if (m_verbose > 0) {
      amrex::Print() << " Doing initial pressure iteration(s) \n";
   }

   for (int lev = 0; lev <= finest_level; ++lev) {
      m_t_old[lev] = m_t_new[lev];
   }

   //----------------------------------------------------------------
   // Initial pressure iterations
   for (int iter = 0; iter < m_init_iter; iter++) {

      if (m_verbose > 0) {
         amrex::Print() << "\n ================   INITIAL ITERATION ["<<iter<<"]   ================ \n";
      }
      int is_init = 1;
      Advance(is_init);

      // Pass new pressure and gp from New to Old
      copyPressNewToOld();

      // Copy back old state
      copyStateOldToNew();
   }
}

void PeleLM::InitFromGridFile(amrex::Real time)
{
  {
     const amrex::BoxArray& ba = MakeBaseGrids();
     DistributionMapping dm(ba);
     MakeNewLevelFromScratch(0, time, ba, dm);
  }
  finest_level = m_initial_ba.size();
  for (int lev = 1; lev <= finest_level; lev++) {
     const amrex::BoxArray ba = m_initial_ba[lev-1];
     DistributionMapping dm(ba);
     MakeNewLevelFromScratch(lev, time, ba, dm);
  }
}

void PeleLM::checkRunParams()
{
#ifdef AMREX_USE_EB
    if (geom[0].IsRZ()) {
        Abort("RZ geometry is not available with EB");
    }
#endif

#if (AMREX_SPACEDIM == 2)
    if (geom[0].IsRZ() && m_phys_bc.lo(0) != 3) {
        Abort("x-low must be 'Symmetry' when using RZ coordinate system");
    }
#endif
}
