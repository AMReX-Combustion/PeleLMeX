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
         amrex::Print() << " with BoxArray " << ba << std::endl;
      }
   }

   // New level factory
#ifdef AMREX_USE_EB
   //TODO
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
                           n_leveldata_new->species, 0);
      fillcoarsepatch_energy(lev, time, n_leveldata_new->rhoh,
                             n_leveldata_new->temp, 0);
      if (m_has_divu) {
         fillcoarsepatch_divu(lev, time, n_leveldata_new->divu,0);
      }
   }

   // Move std::unique_ptr into the PeleLM vector
   m_leveldata_old[lev] = std::move(n_leveldata_old);
   m_leveldata_new[lev] = std::move(n_leveldata_new);
   m_factory[lev]       = std::move(new_fact);

   // DiffusionOp will be recreated
   m_diffusion_op.reset();
   m_diffusionTensor_op.reset();

   // MacProj
#ifdef AMREX_USE_EB
   //TODO
#else
   // TODO: this is weird: at this point AmrCore has not updated finest_level so I need to put a +1
   //       but incflo doesn't need that ...
   macproj.reset(new MacProjector(Geom(0,finest_level+1)));
#endif
}

void PeleLM::RemakeLevel( int lev,
                                amrex::Real time,
                          const amrex::BoxArray& ba,
                          const amrex::DistributionMapping& dm) {
   BL_PROFILE_VAR("PeleLM::RemakeLevel()", RemakeLevel);

   if (m_verbose > 0) {
      Print() << " Remaking level " << lev << "\n";
      if (m_verbose > 2) {
         amrex::Print() << " with BoxArray " << ba << std::endl;
      }
   }

   // New level factory
#ifdef AMREX_USE_EB
   //TODO
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
   fillpatch_velocity(lev, time, n_leveldata_new->velocity, 0);
   fillpatch_gradp(lev, time, n_leveldata_new->gp, 0);
   n_leveldata_new->press.setVal(0.0);

   if (!m_incompressible) {
      fillpatch_density(lev, time, n_leveldata_new->density, 0);
      fillpatch_species(lev, time, n_leveldata_new->species, 0);
      fillpatch_energy(lev, time, n_leveldata_new->rhoh,
                       n_leveldata_new->temp, 0);
      if (m_has_divu) {
         fillpatch_divu(lev, time, n_leveldata_new->divu, 0);
      }
   }

   // Move std::unique_ptr into the PeleLM vector
   m_leveldata_old[lev] = std::move(n_leveldata_old);
   m_leveldata_new[lev] = std::move(n_leveldata_new);
   m_factory[lev]       = std::move(new_fact);

   // DiffusionOp will be recreated
   m_diffusion_op.reset();
   m_diffusionTensor_op.reset();

   // MacProj
#ifdef AMREX_USE_EB
   //TODO
#else
   macproj.reset(new MacProjector(Geom(0,finest_level)));
#endif
}

void PeleLM::ClearLevel(int lev) {
   BL_PROFILE_VAR("PeleLM::ClearLevel()", ClearLevel);

   m_leveldata_old[lev].reset();
   m_leveldata_new[lev].reset();
   m_factory[lev].reset();
   m_diffusion_op.reset();
   m_diffusionTensor_op.reset();
   macproj.reset();
}
