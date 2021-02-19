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
         amrex::Print() << "with BoxArray " << ba << std::endl;
      }
   }

   // Pass Box and Dmap to AmrCore
   SetBoxArray(lev, ba);
   SetDistributionMap(lev, dm);

   // Define the FAB Factory
#ifdef AMREX_USE_EB
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

}

void PeleLM::initData() {
   BL_PROFILE_VAR("PeleLM::initData()", initData);

   if (m_restart_file.empty()) {

      // This is an AmrCore member function which recursively makes new levels
      // with MakeNewLevelFromScratch.
      InitFromScratch(m_cur_time);

      if (m_has_divu) {
         int is_initialization = 1;
         calcDivU(is_initialization,AmrNewTime);
      }

      if (m_do_init_proj) {
         initialProjection();
      }

      initialIterations();

      m_nstep = 0;

      if (m_plot_int > 0 ) {
         WritePlotFile();
      }

   } else {
      // Restart from checkpoint file
   }

}

void PeleLM::initLevelData(int lev) {
   BL_PROFILE_VAR("PeleLM::initLevelData()", initLevelData);

   // Get level data
   auto& ldata = *m_leveldata_new[lev];
   const auto geomdata = geom[lev].data();

   // Pressure and pressure gradients to zero
   ldata.press.setVal(0.0);
   ldata.gp.setVal(0.0);

   // Prob/PMF datas
   ProbParm const* lprobparm = prob_parm.get();
   PmfData const* lpmfdata   = pmf_data_g;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldata.density,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& box = mfi.validbox();
      auto  vel_arr   = ldata.velocity.array(mfi);
      auto  rho_arr   = ldata.density.array(mfi);
      auto  rhoY_arr  = ldata.species.array(mfi);
      auto  rhoH_arr  = ldata.rhoh.array(mfi);
      auto  temp_arr  = ldata.temp.array(mfi);
      auto  aux_arr   = (m_nAux > 0) ? ldata.auxiliaries.array(mfi) : ldata.density.array(mfi);

      amrex::ParallelFor(box,
      [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         pelelm_initdata(i, j, k, vel_arr, rho_arr, rhoY_arr, rhoH_arr, temp_arr, aux_arr,
                         geomdata, *lprobparm, lpmfdata);
      });
   }

   // Initialize thermodynamic pressure
   setThermoPress(lev, AmrNewTime);
}

void PeleLM::initialIterations() {
   BL_PROFILE_VAR("PeleLM::initialIterations()", initialIterations);
}
