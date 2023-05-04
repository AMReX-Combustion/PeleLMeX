#include <PeleLM.H>

using namespace amrex;

void PeleLM::regrid(int lbase,
                    amrex::Real time,
                    bool initial)
{
  BL_PROFILE("PeleLM::regrid()");

  if (!m_doLoadBalance && lbase >= max_level) return;

  if (!m_regrid_file.empty()) {
     regridFromGridFile(lbase, time, initial);

  } else {

     // Load balance base grid
     if (lbase == 0 && m_doLoadBalance && !initial) {

         if ( m_verbose > 0) {
             Print() << " Load balancing level 0 \n";
         }

         int remakeLevel = false;

         computeCosts(0);

         // Use efficiency: average MPI rank cost / max cost
         amrex::Real currentEfficiency = 0.0;
         amrex::Real testEfficiency = 0.0;

         DistributionMapping test_dmap;
         // Build the test dmap, w/o braodcasting
         if (m_loadBalanceMethod == LoadBalanceMethod::SFC) {

             test_dmap = DistributionMapping::makeSFC(*m_costs[0],
                                                      currentEfficiency, testEfficiency,
                                                      false,
                                                      ParallelDescriptor::IOProcessorNumber());

         } else if (m_loadBalanceMethod == LoadBalanceMethod::Knapsack) {

             const amrex::Real navg = static_cast<Real>(grids[0].size()) /
                                      static_cast<Real>(ParallelDescriptor::NProcs());
             const int nmax = static_cast<int>(std::max(std::round(m_loadBalanceKSfactor*navg), std::ceil(navg)));
             test_dmap = DistributionMapping::makeKnapSack(*m_costs[0],
                                                           currentEfficiency, testEfficiency,
                                                           nmax,
                                                           false,
                                                           ParallelDescriptor::IOProcessorNumber());
         }

         // IO proc determine if the test dmap offers significant improvements
         if ((m_loadBalanceEffRatioThreshold > 0.0)
             && (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()))
         {
             remakeLevel = remakeLevel || (testEfficiency > m_loadBalanceEffRatioThreshold*currentEfficiency);
         }
         ParallelDescriptor::Bcast(&remakeLevel, 1, ParallelDescriptor::IOProcessorNumber());

         if ( m_verbose > 1 && remakeLevel) {
             Print() << " Current LoadBalancing efficiency: " << currentEfficiency << "\n"
                     << " Test LoadBalancing efficiency: " << testEfficiency << " \n";
         }

         // Bcast the test dmap and remake level
         if (remakeLevel) {
             Vector<int> pmap;
             if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
                 pmap = test_dmap.ProcessorMap();
             } else {
                 pmap.resize(static_cast<std::size_t>(grids[0].size()));
             }
             ParallelDescriptor::Bcast(pmap.data(), pmap.size(), ParallelDescriptor::IOProcessorNumber());

             if (ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) {
                 test_dmap = DistributionMapping(pmap);
             }

             RemakeLevel(0, time, grids[0], test_dmap);
             SetDistributionMap(0, test_dmap);
         }
     }

     if (lbase >= max_level) return;

     int new_finest;
     Vector<BoxArray> new_grids(finest_level+2);
     MakeNewGrids(lbase, time, new_finest, new_grids);

     BL_ASSERT(new_finest <= finest_level+1);

     bool coarse_ba_changed = false;
     for (int lev = lbase+1; lev <= new_finest; ++lev)
     {
         if (lev <= finest_level) // an old level
         {
             bool ba_changed = (new_grids[lev] != grids[lev]);

             // If initial grid generation or not doing load balancing,
             // remake level only if the grid changed using the default DMap strategy
             if ( initial || !m_doLoadBalance ) {
                 if ( ba_changed || coarse_ba_changed ) {
                     BoxArray level_grids = grids[lev];
                     DistributionMapping level_dmap = dmap[lev];
                     if (ba_changed) {
                         level_grids = new_grids[lev];
                         level_dmap = DistributionMapping(level_grids);
                     }
                     const auto old_num_setdm = num_setdm;
                     RemakeLevel(lev, time, level_grids, level_dmap);
                     SetBoxArray(lev, level_grids);
                     if (old_num_setdm == num_setdm) {
                         SetDistributionMap(lev, level_dmap);
                     }
                 }

             // When doing load balancing, remake level if grid changed or a better dmap
             // can be obtained
             } else {

                 int remakeLevel = false;

                 BoxArray new_ba;
                 DistributionMapping new_dmap;

                 // If the grid changed, let's build a new dmap
                 if (ba_changed) {

                     remakeLevel = true;

                     new_ba = new_grids[lev];
                     new_dmap = DistributionMapping(new_ba);

                     // Get the cost on a LayoutData associated with the new grid
                     LayoutData<Real> new_cost(new_ba,new_dmap);
                     computeCosts(lev, new_cost, m_loadBalanceCost);

                     if (m_loadBalanceMethod == LoadBalanceMethod::SFC) {
                         Vector<Real> costsVec(new_ba.size());
                         ParallelDescriptor::GatherLayoutDataToVector(new_cost, costsVec,
                                                                      ParallelContext::IOProcessorNumberSub());
                         ParallelDescriptor::Bcast(costsVec.data(), costsVec.size(), ParallelContext::IOProcessorNumberSub());
                         Real efficiency;
                         new_dmap = DistributionMapping::makeSFC(costsVec, new_ba, efficiency);

                     } else if (m_loadBalanceMethod == LoadBalanceMethod::Knapsack) {

                         const amrex::Real navg = static_cast<Real>(new_ba.size()) /
                                                  static_cast<Real>(ParallelDescriptor::NProcs());
                         const int nmax = static_cast<int>(std::max(std::round(m_loadBalanceKSfactor*navg), std::ceil(navg)));
                         Vector<Real> costsVec(new_ba.size());
                         ParallelDescriptor::GatherLayoutDataToVector(new_cost, costsVec,
                                                                      ParallelContext::IOProcessorNumberSub());
                         ParallelDescriptor::Bcast(costsVec.data(), costsVec.size(), ParallelContext::IOProcessorNumberSub());
                         Real efficiency;
                         new_dmap = DistributionMapping::makeKnapSack(costsVec, efficiency, nmax);
                     }

                 // Let's see if we can get a better dmap
                 } else {

                     if ( m_verbose > 1) {
                         Print() << " Load balancing level " << lev << "\n";
                     }

                     // Try to build a new dmap with the same old grid
                     new_ba = grids[lev];
                     DistributionMapping test_dmap;

                     computeCosts(lev);

                     // Use efficiency: average MPI rank cost / max cost
                     amrex::Real currentEfficiency = 0.0;
                     amrex::Real testEfficiency = 0.0;

                     // Build the test dmap, w/o braodcasting
                     if (m_loadBalanceMethod == LoadBalanceMethod::SFC) {

                         test_dmap = DistributionMapping::makeSFC(*m_costs[lev],
                                                                  currentEfficiency, testEfficiency,
                                                                  false,
                                                                  ParallelDescriptor::IOProcessorNumber());

                     } else if (m_loadBalanceMethod == LoadBalanceMethod::Knapsack) {

                         const amrex::Real navg = static_cast<Real>(new_ba.size()) /
                                                  static_cast<Real>(ParallelDescriptor::NProcs());
                         const int nmax = static_cast<int>(std::max(std::round(m_loadBalanceKSfactor*navg), std::ceil(navg)));
                         test_dmap = DistributionMapping::makeKnapSack(*m_costs[lev],
                                                                       currentEfficiency, testEfficiency,
                                                                       nmax,
                                                                       false,
                                                                       ParallelDescriptor::IOProcessorNumber());
                     }

                     // IO proc determine if the test dmap offers significant improvements
                     if ((m_loadBalanceEffRatioThreshold > 0.0)
                         && (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()))
                     {
                         remakeLevel = remakeLevel || (testEfficiency > m_loadBalanceEffRatioThreshold*currentEfficiency);
                     }
                     ParallelDescriptor::Bcast(&remakeLevel, 1, ParallelDescriptor::IOProcessorNumber());

                     if ( m_verbose > 1 && remakeLevel) {
                         Print() << " Current LoadBalancing efficiency: " << currentEfficiency << "\n"
                                 << " Test LoadBalancing efficiency: " << testEfficiency << " \n";
                     }

                     // Bcast the test dmap if we plan on remaking the level
                     if (remakeLevel) {
                         Vector<int> pmap;
                         if (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
                             pmap = test_dmap.ProcessorMap();
                         } else {
                             pmap.resize(static_cast<std::size_t>(new_ba.size()));
                         }
                         ParallelDescriptor::Bcast(pmap.data(), pmap.size(), ParallelDescriptor::IOProcessorNumber());

                         if (ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) {
                             test_dmap = DistributionMapping(pmap);
                         }
                         new_dmap = test_dmap;
                     }
                 }

                 if (remakeLevel) {
                     const auto old_num_setdm = num_setdm;
                     RemakeLevel(lev, time, new_ba, new_dmap);
                     SetBoxArray(lev, new_ba);
                     if (old_num_setdm == num_setdm) {
                         SetDistributionMap(lev, new_dmap);
                     }
                 }
             }
             coarse_ba_changed = ba_changed;
         } else {  // a new level, use default DMap strategy
             DistributionMapping new_dmap(new_grids[lev]);
             const auto old_num_setdm = num_setdm;
             MakeNewLevelFromCoarse(lev, time, new_grids[lev], new_dmap);
             SetBoxArray(lev, new_grids[lev]);
             if (old_num_setdm == num_setdm) {
                 SetDistributionMap(lev, new_dmap);
             }
         }
     }

     for (int lev = new_finest+1; lev <= finest_level; ++lev) {
         ClearLevel(lev);
         ClearBoxArray(lev);
         ClearDistributionMap(lev);
     }

     finest_level = new_finest;
  }
}

void PeleLM::MakeNewLevelFromCoarse( int lev,
                                           amrex::Real time,
                                     const amrex::BoxArray& ba,
                                     const amrex::DistributionMapping& dm) {
   BL_PROFILE_VAR("PeleLM::MakeNewLevelFromCoarse()", MakeNewLevelFromCoarse);

   if (m_verbose > 0) {
      Print() << " Making new level " << lev << " from coarse\n";
      if (m_verbose > 2) {
         auto const dx = geom[lev].CellSizeArray();
         Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
         amrex::Print() << " with " << ba.numPts() << " cells, " << ba.size() << " boxes,"
                        << " over " << ba.numPts() * vol / geom[0].ProbSize() * 100 << "% of the domain \n";
      }
      if (m_verbose > 3) {
         amrex::Print() << " with BoxArray " << ba << std::endl;
      }
   }

   // New level factory
#ifdef AMREX_USE_EB
   std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                       {6,6,6},
                                                                       EBSupport::full);
#else
   std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif

   // New leveldatas
   std::unique_ptr<LevelData> n_leveldata_old( new LevelData(ba, dm, *new_fact,
                                                   m_incompressible, m_has_divu,
                                                   m_nAux, m_nGrowState, m_use_soret, m_do_les));

   std::unique_ptr<LevelData> n_leveldata_new( new LevelData(ba, dm, *new_fact,
                                                   m_incompressible, m_has_divu,
                                                   m_nAux, m_nGrowState, m_use_soret, m_do_les));

   // Fill the leveldata_new
   fillcoarsepatch_state(lev, time, n_leveldata_new->state, m_nGrowState);
   fillcoarsepatch_gradp(lev, time, n_leveldata_new->gp, 0);
   n_leveldata_new->press.setVal(0.0);

   if (!m_incompressible) {
      if (m_has_divu) {
         fillcoarsepatch_divu(lev, time, n_leveldata_new->divu,0);
      }
   }

   // Move std::unique_ptr into the PeleLM vector
   m_leveldata_old[lev]  = std::move(n_leveldata_old);
   m_leveldata_new[lev]  = std::move(n_leveldata_new);
   m_factory[lev]        = std::move(new_fact);

   if (m_do_react) {
      std::unique_ptr<LevelDataReact> n_leveldatareact( new LevelDataReact(ba, dm, *m_factory[lev]));
      fillcoarsepatch_reaction(lev, time, n_leveldatareact->I_R, 0);
      n_leveldatareact->functC.setVal(0.0);
      m_leveldatareact[lev] = std::move(n_leveldatareact);
   }

   if (!m_incompressible) {
      // Enforce density / species density consistency
      // only usefull when using cell cons interp
      if (m_regrid_interp_method == 1) {
         setRhoToSumRhoY(lev, AmrNewTime);
      }

      // Recompute temperature from interpolated rhoYs and rhoH
      setTemperature(lev, AmrNewTime);

      // Initialize thermodynamic pressure
      setThermoPress(lev, AmrNewTime);
   }

   if (max_level > 0 && lev != max_level) {
      m_coveredMask[lev].reset(new iMultiFab(ba, dm, 1, 0));
   }
   m_resetCoveredMask = 1;

#ifdef PELE_USE_EFIELD
   m_leveldatanlsolve[lev].reset(new LevelDataNLSolve(ba, dm, *m_factory[lev], m_nGrowState));
   if (m_do_extraEFdiags) {
      m_ionsFluxes[lev].reset(new MultiFab(ba, dm, NUM_IONS*AMREX_SPACEDIM, 0));
   }
   m_precond_op.reset();
#endif

   // Load balance
   m_costs[lev] = std::make_unique<LayoutData<Real>>(ba, dm);

   // DiffusionOp will be recreated
   m_diffusion_op.reset();
   m_mcdiffusion_op.reset();
   m_diffusionTensor_op.reset();

   // Trigger MacProj reset
   m_macProjNeedReset = 1;
   m_extSource[lev].reset(new MultiFab(ba, dm, NVAR, amrex::max(m_nGrowAdv, m_nGrowMAC),
                                       MFInfo(), *m_factory[lev]));
   m_extSource[lev]->setVal(0.);
}

void PeleLM::RemakeLevel( int lev,
                                amrex::Real time,
                          const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm) {
   BL_PROFILE_VAR("PeleLM::RemakeLevel()", RemakeLevel);

   if (m_verbose > 0) {
      Print() << " Remaking level " << lev << "\n";
      if (m_verbose > 2) {
         auto const dx = geom[lev].CellSizeArray();
         Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
         amrex::Print() << " with " << ba.numPts() << " cells," << ba.size() << " boxes,"
                        << " over " << ba.numPts() * vol / geom[0].ProbSize() * 100 << "% of the domain \n";
      }
      if (m_verbose > 3) {
         amrex::Print() << " with BoxArray " << ba << std::endl;
      }
   }

   // New level factory
#ifdef AMREX_USE_EB
   std::unique_ptr<FabFactory<FArrayBox> > new_fact = makeEBFabFactory(geom[lev], ba, dm,
                                                                       {6,6,6},
                                                                       EBSupport::full);
#else
   std::unique_ptr<FabFactory<FArrayBox> > new_fact(new FArrayBoxFactory());
#endif

   // New leveldatas
   std::unique_ptr<LevelData> n_leveldata_old( new LevelData(ba, dm, *new_fact,
                                                   m_incompressible, m_has_divu,
                                                   m_nAux, m_nGrowState, m_use_soret, m_do_les));

   std::unique_ptr<LevelData> n_leveldata_new( new LevelData(ba, dm, *new_fact,
                                                   m_incompressible, m_has_divu,
                                                   m_nAux, m_nGrowState, m_use_soret, m_do_les));

   // Fill the leveldata_new
   fillpatch_state(lev, time, n_leveldata_new->state, m_nGrowState);
   fillpatch_gradp(lev, time, n_leveldata_new->gp, 0);
   n_leveldata_new->press.setVal(0.0);

   if (!m_incompressible) {
      if (m_has_divu) {
         fillpatch_divu(lev, time, n_leveldata_new->divu, 1);
      }
   }

   // Move std::unique_ptr into the PeleLM vector
   m_leveldata_old[lev]  = std::move(n_leveldata_old);
   m_leveldata_new[lev]  = std::move(n_leveldata_new);
   m_factory[lev]        = std::move(new_fact);

   if (m_do_react) {
      std::unique_ptr<LevelDataReact> n_leveldatareact( new LevelDataReact(ba, dm, *m_factory[lev]));
      fillpatch_reaction(lev, time, n_leveldatareact->I_R, 0);
      n_leveldatareact->functC.setVal(0.0);
      m_leveldatareact[lev] = std::move(n_leveldatareact);
   }

   if (max_level > 0 && lev != max_level) {
      m_coveredMask[lev].reset(new iMultiFab(ba, dm, 1, 0));
   }
   m_resetCoveredMask = 1;

   if (!m_incompressible) {
      // Enforce density / species density consistency
      // only usefull when using cell cons interp
      if (m_regrid_interp_method == 1) {
         setRhoToSumRhoY(lev, AmrNewTime);
      }

      // Recompute temperature from interpolated rhoYs and rhoH
      setTemperature(lev, AmrNewTime);

      // Initialize thermodynamic pressure
      setThermoPress(lev, AmrNewTime);
   }

#ifdef PELE_USE_EFIELD
   m_leveldatanlsolve[lev].reset(new LevelDataNLSolve(ba, dm, *m_factory[lev], m_nGrowState));
   if (m_do_extraEFdiags) {
      m_ionsFluxes[lev].reset(new MultiFab(ba, dm, NUM_IONS*AMREX_SPACEDIM, 0));
   }
   m_precond_op.reset();
#endif

   // Load balance
   m_costs[lev] = std::make_unique<LayoutData<Real>>(ba, dm);

   // DiffusionOp will be recreated
   m_diffusion_op.reset();
   m_mcdiffusion_op.reset();
   m_diffusionTensor_op.reset();

   // Trigger MacProj reset
   m_macProjNeedReset = 1;
   m_extSource[lev].reset(new MultiFab(ba, dm, NVAR, amrex::max(m_nGrowAdv, m_nGrowMAC),
                                       MFInfo(), *m_factory[lev]));
   m_extSource[lev]->setVal(0.);
}

void PeleLM::ClearLevel(int lev) {
   BL_PROFILE_VAR("PeleLM::ClearLevel()", ClearLevel);

   m_leveldata_old[lev].reset();
   m_leveldata_new[lev].reset();
   if (m_do_react) m_leveldatareact[lev].reset();
   if (max_level > 0 && lev != max_level) m_coveredMask[lev].reset();
   m_baChem[lev].reset();
   m_dmapChem[lev].reset();
   m_factory[lev].reset();
   m_diffusion_op.reset();
   m_mcdiffusion_op.reset();
   m_diffusionTensor_op.reset();
   macproj.reset();
#ifdef PELE_USE_EFIELD
   m_leveldatanlsolve[lev].reset();
   if (m_do_extraEFdiags) {
      m_ionsFluxes[lev].reset();
   }
#endif
   m_extSource[lev]->clear();

   m_costs[lev].reset();
   m_loadBalanceEff[lev] = -1.0;
}

void PeleLM::computeCosts(int a_lev, LayoutData<Real> &a_costs, int a_costMethod)
{
    if (a_costMethod == LoadBalanceCost::Ncell ) {
        for (MFIter mfi(a_costs, false); mfi.isValid(); ++mfi)
        {
            a_costs[mfi] = mfi.validbox().numPts();
        }
    } else if (a_costMethod == LoadBalanceCost::ChemFunctCallAvg ) {
        MultiFab costMF(a_costs.boxArray(), a_costs.DistributionMap(), 1, 0);
        fillpatch_chemFunctCall(a_lev, m_cur_time, costMF, 0);
        for (MFIter mfi(costMF, false); mfi.isValid(); ++mfi)
        {
            a_costs[mfi] = costMF[mfi].sum<RunOn::Device>(mfi.validbox(),0) / mfi.validbox().numPts();
        }
    } else if (a_costMethod == LoadBalanceCost::ChemFunctCallMax ) {
        MultiFab costMF(a_costs.boxArray(), a_costs.DistributionMap(), 1, 0);
        fillpatch_chemFunctCall(a_lev, m_cur_time, costMF, 0);
        for (MFIter mfi(costMF, false); mfi.isValid(); ++mfi)
        {
            a_costs[mfi] = costMF[mfi].max<RunOn::Device>(mfi.validbox(),0);
        }
    } else if (a_costMethod == LoadBalanceCost::ChemFunctCallSum ) {
        MultiFab costMF(a_costs.boxArray(), a_costs.DistributionMap(), 1, 0);
        fillpatch_chemFunctCall(a_lev, m_cur_time, costMF, 0);
        for (MFIter mfi(costMF, false); mfi.isValid(); ++mfi)
        {
            a_costs[mfi] = costMF[mfi].sum<RunOn::Device>(mfi.validbox(),0);
        }
    } else if (a_costMethod == LoadBalanceCost::UserDefinedDerivedAvg) {
        MultiFab costMF(a_costs.boxArray(), a_costs.DistributionMap(), 1, 0);
        costMF.setVal(0.0);
        std::unique_ptr<MultiFab> mf;
        mf = derive("derUserDefined", m_cur_time, a_lev, 0);
        costMF.ParallelCopy(*mf,0,0,1);
        for (MFIter mfi(costMF, false); mfi.isValid(); ++mfi)
        {
            a_costs[mfi] = costMF[mfi].sum<RunOn::Device>(mfi.validbox(),0) / mfi.validbox().numPts();
        }
    } else if (a_costMethod == LoadBalanceCost::UserDefinedDerivedSum) {
        MultiFab costMF(a_costs.boxArray(), a_costs.DistributionMap(), 1, 0);
        costMF.setVal(0.0);
        std::unique_ptr<MultiFab> mf;
        mf = derive("derUserDefined", m_cur_time, a_lev, 0);
        costMF.ParallelCopy(*mf,0,0,1);
        for (MFIter mfi(costMF, false); mfi.isValid(); ++mfi)
        {
            a_costs[mfi] = costMF[mfi].sum<RunOn::Device>(mfi.validbox(),0);
        }
    } else {
        Abort(" Unknown cost estimate method !");
    }
}

void PeleLM::computeCosts(int a_lev)
{
    computeCosts(a_lev, *m_costs[a_lev], m_loadBalanceCost);
}

void PeleLM::resetMacProjector()
{
   // If nothing has changed, just go back
   if ( !m_macProjNeedReset && ( m_macProjOldSize == finest_level+1 ) ) {
      return;
   }

   // MacProj
#ifdef AMREX_USE_EB
   macproj.reset(new Hydro::MacProjector(Geom(0,finest_level),
                                         MLMG::Location::FaceCentroid,  // Location of mac velocity
                                         MLMG::Location::FaceCentroid,  // Location of beta
                                         MLMG::Location::CellCenter) ); // Location of solution variable phi
#else
   macproj.reset(new Hydro::MacProjector(Geom(0,finest_level)));
#endif

   // Store the old MacProj size and switch off reset flag
   m_macProjOldSize = finest_level+1;
   m_macProjNeedReset = 0;
}

void PeleLM::regridFromGridFile(int lbase,
                                amrex::Real time,
                                bool /*initial*/)
{
  int new_finest = m_regrid_ba.size();
  Vector<BoxArray> new_grids(finest_level+2);
  BL_ASSERT(new_finest <= finest_level+1);

  bool coarse_ba_changed = false;
  for (int lev = lbase+1; lev <= new_finest; ++lev) {
     new_grids[lev] = m_regrid_ba[lev-1];
     if (lev <= finest_level) { // an old level
        bool ba_changed = (new_grids[lev] != grids[lev]);
        if (ba_changed || coarse_ba_changed) {
           BoxArray level_grids = grids[lev];
           DistributionMapping level_dmap = dmap[lev];
           if (ba_changed) {
              level_grids = new_grids[lev];
              level_dmap = DistributionMapping(level_grids);
           }
           const auto old_num_setdm = num_setdm;
           RemakeLevel(lev, time, level_grids, level_dmap);
           SetBoxArray(lev, level_grids);
           if (old_num_setdm == num_setdm) {
              SetDistributionMap(lev, level_dmap);
           }
        }
        coarse_ba_changed = ba_changed;;
      } else { // a new level
        DistributionMapping new_dmap(new_grids[lev]);
        const auto old_num_setdm = num_setdm;
        MakeNewLevelFromCoarse(lev, time, new_grids[lev], new_dmap);
        SetBoxArray(lev, new_grids[lev]);
        if (old_num_setdm == num_setdm) {
           SetDistributionMap(lev, new_dmap);
        }
     }
  }
  for (int lev = new_finest+1; lev <= finest_level; ++lev) {
     ClearLevel(lev);
     ClearBoxArray(lev);
     ClearDistributionMap(lev);
  }

  finest_level = new_finest;
}
