#include <PeleLM.H>
#include <pmf_data.H>
#include <pelelm_prob.H>

using namespace amrex;

void PeleLM::Init() {
   BL_PROFILE_VAR("PeleLM::Init()", Init);

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
      if (m_verbose > 2) {
         amrex::Print() << " with BoxArray " << ba << std::endl;
      }
   }

   // Pass Box and Dmap to AmrCore
   SetBoxArray(lev, ba);
   SetDistributionMap(lev, dm);

   // Define the FAB Factory
#ifdef AMREX_USE_EB
   m_factory[lev] = makeEBFabFactory(geom[lev], grids[lev], dmap[lev],
                                     {nghost_eb_basic(),
                                      nghost_eb_volume(),
                                      nghost_eb_full()},
                                     EBSupport::full);
#else
   m_factory[lev].reset(new FArrayBoxFactory());
#endif

   // Initialize the LevelData
   m_leveldata_old[lev].reset(new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                            m_incompressible, m_has_divu,
                                            m_nAux, m_nGrowState, m_nGrowMAC));
   m_leveldata_new[lev].reset(new LevelData(grids[lev], dmap[lev], *m_factory[lev],
                                            m_incompressible, m_has_divu,
                                            m_nAux, m_nGrowState, m_nGrowMAC));

   if (max_level > 0 && lev != max_level) {
      m_coveredMask[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, 0));
      m_resetCoveredMask = 1;
   }
   if (m_do_react) {
      m_leveldatareact[lev].reset(new LevelDataReact(grids[lev], dmap[lev], *m_factory[lev]));
      m_leveldatareact[lev]->functC.setVal(0.0);
   }

#ifdef PLM_USE_EFIELD
   int nGrowNL = 1;
   m_leveldatanlsolve[lev].reset(new LevelDataNLSolve(grids[lev], dmap[lev], *m_factory[lev], 1));
#endif

   // Fill the initial solution (if not restarting)
   if (m_restart_file.empty()) {
      initLevelData(lev);
   }

   // Times
   m_t_new[lev] = time;
   m_t_old[lev] = time - 1.0e200;

   // Mac projector
#if AMREX_USE_EB
   macproj.reset(new MacProjector(Geom(0,finest_level),
                                  MLMG::Location::FaceCentroid,  // Location of mac velocity
                                  MLMG::Location::FaceCentroid,  // Location of beta
                                  MLMG::Location::CellCenter) ); // Location of solution variable phi
#else
   macproj.reset(new MacProjector(Geom(0,finest_level)));
#endif
   m_macProjOldSize = finest_level+1;
}

void PeleLM::initData() {
   BL_PROFILE_VAR("PeleLM::initData()", initData);

   if (m_restart_file.empty()) {

      //----------------------------------------------------------------
      // This is an AmrCore member function which recursively makes new levels
      // with MakeNewLevelFromScratch.
      InitFromScratch(m_cur_time);
      resetCoveredMask();

      //----------------------------------------------------------------
      // AverageDown and FillPatch the NewState
      averageDownState(AmrNewTime);
      fillPatchState(AmrNewTime);

#ifdef PLM_USE_EFIELD
      poissonSolveEF(AmrNewTime);
#endif

      // Post data Init time step estimate
      // TODO : this estimate is probably useless
      int is_init = 1;
      Real dtInit = computeDt(is_init,AmrNewTime);
      Print() << " Initial dt: " << dtInit << "\n";

      // Subcycling IAMR/PeleLM first does a projection without divU
      // whch makes the dt for evaluating I_R better
      int has_divu_save = m_has_divu;
      m_has_divu = 0;
      initialProjection();
      m_has_divu = has_divu_save;

      // Post data Init time step estimate
      dtInit = computeDt(is_init,AmrNewTime);
      Print() << " Initial dt: " << dtInit << "\n";

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
                  MultiFab Forcing(grids[lev],dmap[lev],NUM_SPECIES+1,0);
                  Forcing.setVal(0.0);

                  if (lev != finest_level) {
                     // On all but the finest level, average down I_R
                     // and use advanceChemistry with chem BoxArray
                     std::unique_ptr<MultiFab> avgDownIR;
                     avgDownIR.reset(new MultiFab(grids[lev],dmap[lev],NUM_SPECIES,0));
                     avgDownIR->setVal(0.0);
                     auto ldataRFine_p   = getLevelDataReactPtr(lev+1);
#ifdef AMREX_USE_EB
                     EB_average_down(ldataRFine_p->I_R,
                                     *avgDownIR,
                                     0,NUM_SPECIES,refRatio(lev));
#else
                     average_down(ldataRFine_p->I_R,
                                  *avgDownIR,
                                  0,NUM_SPECIES,refRatio(lev));
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
            // TODO: closed_chamber correction

            initialProjection();
         }
      }

      //----------------------------------------------------------------
      // Initial pressure iterations
      initialIterations();

      m_nstep = 0;

      if (m_plot_int > 0 ) {
         WritePlotFile();
      }
      if (m_check_int > 0 ) {
         WriteCheckPointFile();
      }

   } else {
      //----------------------------------------------------------------
      // Read starting configuration from chk file.
      ReadCheckPointFile();

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
   ProbParm const* lprobparm = prob_parm.get();
   PmfData const* lpmfdata   = pmf_data_g;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldata_p->velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      FArrayBox DummyFab(bx,1);
      auto  const &vel_arr   = ldata_p->velocity.array(mfi);
      auto  const &rho_arr   = (m_incompressible) ? DummyFab.array() : ldata_p->density.array(mfi);
      auto  const &rhoY_arr  = (m_incompressible) ? DummyFab.array() : ldata_p->species.array(mfi);
      auto  const &rhoH_arr  = (m_incompressible) ? DummyFab.array() : ldata_p->rhoh.array(mfi);
      auto  const &temp_arr  = (m_incompressible) ? DummyFab.array() : ldata_p->temp.array(mfi);
      auto  const &aux_arr   = (m_nAux > 0) ? ldata_p->auxiliaries.array(mfi) : DummyFab.array();
#ifdef PLM_USE_EFIELD
      auto  const &ne_arr    = ldata_p->nE.array(mfi);
      auto  const &phiV_arr  = ldata_p->phiV.array(mfi);
#endif
      amrex::ParallelFor(bx, [=,m_incompressible=m_incompressible]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         pelelm_initdata(i, j, k, m_incompressible, vel_arr, rho_arr,
                         rhoY_arr, rhoH_arr, temp_arr, aux_arr,
#ifdef PLM_USE_EFIELD
                         ne_arr, phiV_arr,  
#endif
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
