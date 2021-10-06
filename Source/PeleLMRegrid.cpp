#include <PeleLM.H>

using namespace amrex;

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
         amrex::Print() << " with " << ba.numPts() << " cells,"
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
                                                   m_nAux, m_nGrowState, m_nGrowMAC));

   std::unique_ptr<LevelData> n_leveldata_new( new LevelData(ba, dm, *new_fact,
                                                   m_incompressible, m_has_divu,
                                                   m_nAux, m_nGrowState, m_nGrowMAC));

   // Fill the leveldata_new
   fillcoarsepatch_velocity(lev, time, n_leveldata_new->velocity, 0);
   fillcoarsepatch_gradp(lev, time, n_leveldata_new->gp, 0);
   n_leveldata_new->press.setVal(0.0);

   if (!m_incompressible) {
      fillcoarsepatch_mass(lev, time, n_leveldata_new->density,
                           n_leveldata_new->species, 1);
      fillcoarsepatch_energy(lev, time, n_leveldata_new->rhoh,
                             n_leveldata_new->temp, 1);
      if (m_has_divu) {
         fillcoarsepatch_divu(lev, time, n_leveldata_new->divu,0);
      }
#ifdef PLM_USE_EFIELD
      fillcoarsepatch_phiV(lev, time, n_leveldata_new->phiV,0);
      fillcoarsepatch_nE(lev, time, n_leveldata_new->nE,0);
#endif
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
      // Initialize thermodynamic pressure
      setThermoPress(lev, AmrNewTime);
   }

   if (max_level > 0 && lev != max_level) {
      m_coveredMask[lev].reset(new iMultiFab(ba, dm, 1, 0));
   }
   m_resetCoveredMask = 1;

#ifdef PLM_USE_EFIELD
   m_leveldatanlsolve[lev].reset(new LevelDataNLSolve(ba, dm, *m_factory[lev], 1));
   m_precond_op.reset();
#endif

   // DiffusionOp will be recreated
   m_diffusion_op.reset();
   m_diffusionTensor_op.reset();

   // Trigger MacProj reset
   m_macProjNeedReset = 1;
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
         amrex::Print() << " with " << ba.numPts() << " cells,"
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
                                                   m_nAux, m_nGrowState, m_nGrowMAC));

   std::unique_ptr<LevelData> n_leveldata_new( new LevelData(ba, dm, *new_fact,
                                                   m_incompressible, m_has_divu,
                                                   m_nAux, m_nGrowState, m_nGrowMAC));

   // Fill the leveldata_new
   fillpatch_velocity(lev, time, n_leveldata_new->velocity, m_nGrowState);
   fillpatch_gradp(lev, time, n_leveldata_new->gp, 0);
   n_leveldata_new->press.setVal(0.0);

   if (!m_incompressible) {
      fillpatch_density(lev, time, n_leveldata_new->density, m_nGrowState);
      fillpatch_species(lev, time, n_leveldata_new->species, m_nGrowState);
      fillpatch_energy(lev, time, n_leveldata_new->rhoh,
                       n_leveldata_new->temp, m_nGrowState);
      if (m_has_divu) {
         fillpatch_divu(lev, time, n_leveldata_new->divu, 1);
      }
#ifdef PLM_USE_EFIELD
      fillpatch_phiV(lev, time, n_leveldata_new->phiV,m_nGrowState);
      fillpatch_nE(lev, time, n_leveldata_new->nE,m_nGrowState);
#endif
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
      // Initialize thermodynamic pressure
      setThermoPress(lev, AmrNewTime);
   }

#ifdef PLM_USE_EFIELD
   m_leveldatanlsolve[lev].reset(new LevelDataNLSolve(ba, dm, *m_factory[lev], 1));
   m_precond_op.reset();
#endif

   // DiffusionOp will be recreated
   m_diffusion_op.reset();
   m_diffusionTensor_op.reset();

   // Trigger MacProj reset
   m_macProjNeedReset = 1;
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
   m_diffusionTensor_op.reset();
   macproj.reset();
#ifdef PLM_USE_EFIELD
   m_leveldatanlsolve[lev].reset();
#endif
}

void PeleLM::resetMacProjector()
{
   // If nothing has changed, just go back
   if ( !m_macProjNeedReset && ( m_macProjOldSize == finest_level+1 ) ) {
      return;
   }

   // MacProj
#ifdef AMREX_USE_EB
   macproj.reset(new MacProjector(Geom(0,finest_level),
                                  MLMG::Location::FaceCentroid,  // Location of mac velocity
                                  MLMG::Location::FaceCentroid,  // Location of beta
                                  MLMG::Location::CellCenter) ); // Location of solution variable phi
#else
   macproj.reset(new MacProjector(Geom(0,finest_level)));
#endif

   // Store the old MacProj size and switch off reset flag
   m_macProjOldSize = finest_level+1;
   m_macProjNeedReset = 0;
}
