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
   //TODO
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

   // Fill the initial solution (if not restarting)
   if (m_restart_file.empty()) {
      initLevelData(lev);
   }

   // Times
   m_t_new[lev] = time;
   m_t_old[lev] = time - 1.0e200;

   // Mac projector
#if AMREX_USE_EB
#else
   macproj.reset(new MacProjector(Geom(0,finest_level)));
#endif

}

void PeleLM::initData() {
   BL_PROFILE_VAR("PeleLM::initData()", initData);

   if (m_restart_file.empty()) {

      // This is an AmrCore member function which recursively makes new levels
      // with MakeNewLevelFromScratch.
      InitFromScratch(m_cur_time);

      // FillPatch the NewState
      fillPatchState(AmrNewTime);

      // Initial velocity projection iterations
      if (m_do_init_proj) {
         for (int iter = 0; iter < m_numDivuIter; iter++) {
            if (m_has_divu) {
               int is_initialization = 1;
               int do_avgDown = 1;
               calcDivU(is_initialization,do_avgDown,AmrNewTime);
            }
            // TODO: closed_chamber correction

            initialProjection();
         }
      }

      // Initial pressure iterations
      initialIterations();

      m_nstep = 0;

      if (m_plot_int > 0 ) {
         WritePlotFile();
      }

   } else {
      // TODO Restart from checkpoint file
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
      amrex::ParallelFor(bx, [=,m_incompressible=m_incompressible]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         pelelm_initdata(i, j, k, m_incompressible, vel_arr, rho_arr, 
                         rhoY_arr, rhoH_arr, temp_arr, aux_arr,
                         geomdata, *lprobparm, lpmfdata);
      });
   }

   if (!m_incompressible) {
      // Initialize thermodynamic pressure
      setThermoPress(lev, AmrNewTime);
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

   // Initial pressure iterations
   for (int iter = 0; iter < m_init_iter; iter++) {

      if (m_verbose > 0) {
         amrex::Print() << "\n ================   INITIAL INTERATION ["<<iter<<"]   ================ \n";
      }
      int is_init = 1;
      Advance(is_init);

      // Pass new pressure and gp from New to Old
      copyPressNewToOld();

      // Copy back old state
      copyStateOldToNew();
   }
}
