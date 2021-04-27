#include <PeleLM.H>
#include <pmf_data.H>
#include <PeleLMBCfill.H>
#include <AMReX_FillPatchUtil.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

// Conversion from physBC to fieldBC maps
// Components are  Interior, Inflow, Outflow, Symmetry, &
// SlipWallAdiab, NoSlipWallAdiab, SlipWallIsoTherm, NoSlipWallIsoTherm.

int
norm_vel_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD
};

int
tang_vel_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_ODD, REFLECT_EVEN, REFLECT_ODD
};

int
density_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

int
species_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

int
rhoh_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

int
temp_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

int
divu_bc[] =
{
  INT_DIR, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

// Following incflo rather than IAMR here
int
force_bc[] =
{
  BCType::int_dir, BCType::foextrap, BCType::foextrap , BCType::foextrap, BCType::foextrap, BCType::foextrap,
  BCType::foextrap, BCType::foextrap
};

#ifdef PLM_USE_EFIELD
int
nE_bc[] =
{
  INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN, EXT_DIR, EXT_DIR
};

int
phiV_bc[] =
{
        INT_DIR, EXT_DIR, REFLECT_EVEN
};
#endif

void PeleLM::setBoundaryConditions() {

   // Initialize the BCRecs
   m_bcrec_state.resize(NVAR);
   int sizeForceBC = std::max(AMREX_SPACEDIM,NUM_SPECIES+2);
   m_bcrec_force.resize(sizeForceBC);

   // Convert m_phys_bc into field BCs
   // Get m_phys_bc
   const int* lo_bc = m_phys_bc.lo();
   const int* hi_bc = m_phys_bc.hi();

   // Velocity
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      for (int idim2 = 0; idim2 < AMREX_SPACEDIM; idim2++) {
         if ( idim == idim2 ) {
            m_bcrec_state[VELX+idim].setLo(idim2,norm_vel_bc[lo_bc[idim2]]);
            m_bcrec_state[VELX+idim].setHi(idim2,norm_vel_bc[hi_bc[idim2]]);
         } else {
            m_bcrec_state[VELX+idim].setLo(idim2,tang_vel_bc[lo_bc[idim2]]);
            m_bcrec_state[VELX+idim].setHi(idim2,tang_vel_bc[hi_bc[idim2]]);
         }
      }
   }  

   // General forces: use int_dir in interior and foextrap otherwise
   for (int i = 0; i < sizeForceBC; i++) {
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_force[i].setLo(idim,force_bc[lo_bc[idim]]);
         m_bcrec_force[i].setHi(idim,force_bc[hi_bc[idim]]);
      }
   }

   if (!m_incompressible) {
      // Density
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_state[DENSITY].setLo(idim,density_bc[lo_bc[idim]]);
         m_bcrec_state[DENSITY].setHi(idim,density_bc[hi_bc[idim]]);
      }

      // Species
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         for ( int n = 0; n < NUM_SPECIES; n++) {
            m_bcrec_state[FIRSTSPEC+n].setLo(idim,density_bc[lo_bc[idim]]);
            m_bcrec_state[FIRSTSPEC+n].setHi(idim,density_bc[hi_bc[idim]]);
         }
      }

      // Enthalpy
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_state[RHOH].setLo(idim,rhoh_bc[lo_bc[idim]]);
         m_bcrec_state[RHOH].setHi(idim,rhoh_bc[hi_bc[idim]]);
      }

      // Temperature
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_state[TEMP].setLo(idim,temp_bc[lo_bc[idim]]);
         m_bcrec_state[TEMP].setHi(idim,temp_bc[hi_bc[idim]]);
      }

      // rhoRT: reflect even on all but interior bndy
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_state[RHORT].setLo(idim,divu_bc[lo_bc[idim]]);
         m_bcrec_state[RHORT].setHi(idim,divu_bc[hi_bc[idim]]);
      }

      // divU
      if (m_has_divu) {
         for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            m_bcrec_divu.setLo(idim,divu_bc[lo_bc[idim]]);
            m_bcrec_divu.setHi(idim,divu_bc[hi_bc[idim]]);
         }
      }

#ifdef PLM_USE_EFIELD
      // nE
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_state[NE].setLo(idim,nE_bc[lo_bc[idim]]);
         m_bcrec_state[NE].setHi(idim,nE_bc[hi_bc[idim]]);
      }

      // Get m_phiV_bc
      const int* lo_phibc = m_phiV_bc.lo();
      const int* hi_phibc = m_phiV_bc.hi();
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_bcrec_state[PHIV].setLo(idim,phiV_bc[lo_phibc[idim]]);
         m_bcrec_state[PHIV].setHi(idim,phiV_bc[hi_phibc[idim]]);
      }

      // Hack charged species BCs
      int FIRSTION = FIRSTSPEC + NUM_SPECIES - NUM_IONS;
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         for ( int n = 0; n < NUM_IONS; n++) {
            auto const bcIonSave = m_bcrec_state[FIRSTION+n];
            m_bcrec_state[FIRSTION+n] = hackBCChargedParticle(zk[FIRSTION+n], bcIonSave);
         }
      }
      // Need to hack nE too actually ...
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         auto const bcnESave = m_bcrec_state[NE];
         m_bcrec_state[NE] = hackBCChargedParticle(-1.0, bcnESave);
      }
#endif
   }

   // Aux
   //TODO
}

Vector<BCRec>
PeleLM::fetchBCRecArray(int scomp, int ncomp) {
   Vector<BCRec> bc(ncomp);
   for (int comp = 0; comp < ncomp; comp++) {
      bc[comp] = m_bcrec_state[scomp+comp];
   }
   return bc;
}


//-----------------------------------------------------------------------------
// The following work directly on the leveldata

// Fill the entire class state at once
void PeleLM::fillPatchState(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchState()", fillPatchState);
   for (int lev = 0; lev <= finest_level; lev++) {
      fillPatchState(lev,a_time);
   }
}

// Fill the a given level class state
void PeleLM::fillPatchState(int lev, const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchStateLev()", fillPatchStateLev);

   auto ldata_p = getLevelDataPtr(lev,a_time);
   Real time = getTime(lev, a_time);

   fillpatch_velocity(lev, time, ldata_p->velocity, m_nGrowState);
   if (!m_incompressible) {
      fillpatch_density(lev, time, ldata_p->density, m_nGrowState);
      fillpatch_species(lev, time, ldata_p->species, m_nGrowState);
      fillpatch_energy(lev, time, ldata_p->rhoh, ldata_p->temp, m_nGrowState);
      fillpatch_thermoPress(lev, time, ldata_p->rhoRT, m_nGrowState);
      if (m_has_divu) {
         fillpatch_divu(lev, time, ldata_p->divu, ldata_p->divu.nGrow());
      }
#ifdef PLM_USE_EFIELD
      fillpatch_nE(lev, time, ldata_p->nE, m_nGrowState);
      fillpatch_phiV(lev, time, ldata_p->phiV, m_nGrowState);
#endif
   }
   //TODO Aux
}

// Fill a state components
void PeleLM::fillPatchDensity(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchDensity()", fillPatchDensity);
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,a_time);
      Real time = getTime(lev, a_time);
      fillpatch_density(lev, time, ldata_p->density, m_nGrowState);
   }
}

void PeleLM::fillPatchSpecies(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchSpecies()", fillPatchSpecies);
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,a_time);
      Real time = getTime(lev, a_time);
      fillpatch_species(lev, time, ldata_p->species, m_nGrowState);
   }
}

void PeleLM::fillPatchTemp(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchTemp()", fillPatchTemp);
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,a_time);
      Real time = getTime(lev, a_time);
      fillpatch_energy(lev, time, ldata_p->rhoh, ldata_p->temp, m_nGrowState);
   }
}

void PeleLM::fillPatchPhiV(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchPhiV()", fillPatchPhiV);
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,a_time);
      Real time = getTime(lev, a_time);
      fillpatch_phiV(lev, time, ldata_p->phiV, m_nGrowState);
   }
}
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// The following return a fillpatched MF ptr at a given level
// Fill the entire state at once
std::unique_ptr<MultiFab>
PeleLM::fillPatchState(int lev, Real a_time, int nGrow) {
   BL_PROFILE_VAR("PeleLM::fillPatchState()", fillPatchState);

   std::unique_ptr<MultiFab> mf;
   mf.reset(new MultiFab(grids[lev], dmap[lev], NVAR, nGrow));

   MultiFab vel(*mf, amrex::make_alias, VELX, AMREX_SPACEDIM);
   fillpatch_velocity(lev, a_time, vel, nGrow);
   if (!m_incompressible) {
      MultiFab rho(*mf, amrex::make_alias, DENSITY, 1);
      fillpatch_density(lev, a_time, rho, nGrow);
      MultiFab rhoY(*mf, amrex::make_alias, FIRSTSPEC, NUM_SPECIES);
      fillpatch_species(lev, a_time, rhoY, nGrow);
      MultiFab rhoh(*mf, amrex::make_alias, RHOH, 1);
      MultiFab temp(*mf, amrex::make_alias, TEMP, 1);
      fillpatch_energy(lev, a_time, rhoh, temp, nGrow);
      MultiFab rhoRT(*mf, amrex::make_alias, RHORT, 1);
      fillpatch_thermoPress(lev, a_time, rhoRT, nGrow);
#ifdef PLM_USE_EFIELD
      MultiFab nE(*mf, amrex::make_alias, NE, 1);
      fillpatch_nE(lev, a_time, nE, nGrow);
      MultiFab phiV(*mf, amrex::make_alias, PHIV, 1);
      fillpatch_phiV(lev, a_time, phiV, nGrow);
#endif
   }
   //TODO Aux

   return mf;
}
//-----------------------------------------------------------------------------

// Fill the velocity
void PeleLM::fillpatch_velocity(int lev,
                                const amrex::Real a_time,
                                amrex::MultiFab &a_vel,
                                int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> bndry_func(geom[lev], fetchBCRecArray(VELX,AMREX_SPACEDIM),
                                                                     PeleLMCCFillExtDirVel{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_vel, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->velocity),&(m_leveldata_new[lev]->velocity)},
                           {m_t_old[lev], m_t_new[lev]},0,0,AMREX_SPACEDIM,geom[lev], bndry_func, 0);
   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> crse_bndry_func(geom[lev-1], fetchBCRecArray(VELX,AMREX_SPACEDIM),
                                                                          PeleLMCCFillExtDirVel{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> fine_bndry_func(geom[lev], fetchBCRecArray(VELX,AMREX_SPACEDIM),
                                                                          PeleLMCCFillExtDirVel{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_vel, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->velocity),&(m_leveldata_new[lev-1]->velocity)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->velocity),&(m_leveldata_new[lev]->velocity)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(VELX,AMREX_SPACEDIM), 0);
   }
}

// Fill the density
void PeleLM::fillpatch_density(int lev,
                               const amrex::Real a_time,
                               amrex::MultiFab &a_density,
                               int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDens>> bndry_func_rho(geom[lev], fetchBCRecArray(DENSITY,1),
                                                                          PeleLMCCFillExtDirDens{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_density, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->density),&(m_leveldata_new[lev]->density)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func_rho, 0);

   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDens>> crse_bndry_func_rho(geom[lev-1], fetchBCRecArray(DENSITY,1), 
                                                                               PeleLMCCFillExtDirDens{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDens>> fine_bndry_func_rho(geom[lev], fetchBCRecArray(DENSITY,1),
                                                                               PeleLMCCFillExtDirDens{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_density, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->density),&(m_leveldata_new[lev-1]->density)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->density),&(m_leveldata_new[lev]->density)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func_rho,0,fine_bndry_func_rho,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(DENSITY,1), 0);
   }
}

// Fill the mass fractions
void PeleLM::fillpatch_species(int lev,
                               const amrex::Real a_time,
                               amrex::MultiFab &a_species,
                               int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {

      // Species
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirSpec>> bndry_func(geom[lev], fetchBCRecArray(FIRSTSPEC,NUM_SPECIES),
                                                                      PeleLMCCFillExtDirSpec{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_species, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->species),&(m_leveldata_new[lev]->species)},
                           {m_t_old[lev], m_t_new[lev]},0,0,NUM_SPECIES,geom[lev], bndry_func, 0);
   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      // Species
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirSpec>> crse_bndry_func(geom[lev-1], fetchBCRecArray(FIRSTSPEC,NUM_SPECIES),
                                                                           PeleLMCCFillExtDirSpec{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirSpec>> fine_bndry_func(geom[lev], fetchBCRecArray(FIRSTSPEC,NUM_SPECIES),
                                                                           PeleLMCCFillExtDirSpec{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_species, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->species),&(m_leveldata_new[lev-1]->species)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->species),&(m_leveldata_new[lev]->species)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, NUM_SPECIES, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(FIRSTSPEC,NUM_SPECIES), 0);
   }
}

// Fill rhoH and temperature
void PeleLM::fillpatch_energy(int lev,
                              const amrex::Real a_time,
                              amrex::MultiFab &a_rhoh,
                              amrex::MultiFab &a_temp,
                              int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {

      // rhoH
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirRhoH>> bndry_func_rhoh(geom[lev], fetchBCRecArray(RHOH,1),
                                                                          PeleLMCCFillExtDirRhoH{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_rhoh, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->rhoh),&(m_leveldata_new[lev]->rhoh)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func_rhoh, 0);

      // temperature
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirTemp>> bndry_func(geom[lev], fetchBCRecArray(TEMP,1),
                                                                      PeleLMCCFillExtDirTemp{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_temp, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->temp),&(m_leveldata_new[lev]->temp)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func, 0);
   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      // rhoH
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirRhoH>> crse_bndry_func_rhoh(geom[lev-1], fetchBCRecArray(RHOH,1),
                                                                                PeleLMCCFillExtDirRhoH{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirRhoH>> fine_bndry_func_rhoh(geom[lev], fetchBCRecArray(RHOH,1),
                                                                                PeleLMCCFillExtDirRhoH{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_rhoh, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->rhoh),&(m_leveldata_new[lev-1]->rhoh)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->rhoh),&(m_leveldata_new[lev]->rhoh)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func_rhoh,0,fine_bndry_func_rhoh,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(RHOH,1), 0);

      // temperature
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirTemp>> crse_bndry_func(geom[lev-1], fetchBCRecArray(TEMP,1),
                                                                           PeleLMCCFillExtDirTemp{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirTemp>> fine_bndry_func(geom[lev], fetchBCRecArray(TEMP,1),
                                                                           PeleLMCCFillExtDirTemp{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_temp, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->temp),&(m_leveldata_new[lev-1]->temp)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->temp),&(m_leveldata_new[lev]->temp)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(TEMP,1), 0);

   }
}

// Fill the thermodynamic pressure
void PeleLM::fillpatch_thermoPress(int lev,
                                   const amrex::Real a_time,
                                   amrex::MultiFab &a_rhoRT,
                                   int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> bndry_func_rhoRT(geom[lev], fetchBCRecArray(RHORT,1),
                                                                             PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchSingleLevel(a_rhoRT, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->rhoRT),&(m_leveldata_new[lev]->rhoRT)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func_rhoRT, 0);

   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func_rhoRT(geom[lev-1], fetchBCRecArray(RHORT,1), 
                                                                                  PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func_rhoRT(geom[lev], fetchBCRecArray(RHORT,1),
                                                                                  PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchTwoLevels(a_rhoRT, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->rhoRT),&(m_leveldata_new[lev-1]->rhoRT)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->rhoRT),&(m_leveldata_new[lev]->rhoRT)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func_rhoRT,0,fine_bndry_func_rhoRT,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(DENSITY,1), 0);
   }
}

// Fill the divU
void PeleLM::fillpatch_divu(int lev,
                            const amrex::Real a_time,
                            amrex::MultiFab &a_divu,
                            int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> bndry_func(geom[lev], {m_bcrec_divu},
                                                                     PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchSingleLevel(a_divu, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->divu),&(m_leveldata_new[lev]->divu)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func, 0);
   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func(geom[lev-1], {m_bcrec_divu},
                                                                            PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func(geom[lev], {m_bcrec_divu},
                                                                            PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchTwoLevels(a_divu, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->divu),&(m_leveldata_new[lev-1]->divu)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->divu),&(m_leveldata_new[lev]->divu)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_divu}, 0);
   }
}

// Fillpatch a vector of forces:
// -> actually only modifies the ghost cells : fillBoundary, C/F interp, foextrap on domain BCs
void PeleLM::fillpatch_forces(Real a_time,
                              Vector<MultiFab*> const &a_force,
                              int nGrowForce)
{
   AMREX_ASSERT(a_force[0]->nComp() <= m_bcrec_force.size());
   const int nComp = a_force[0]->nComp();
   ProbParm const* lprobparm = prob_parm.get();

   int lev = 0;
   {
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy> > bndry_func(geom[lev], {m_bcrec_force},
                                                                        PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchSingleLevel(*a_force[lev],IntVect(nGrowForce),a_time,{a_force[lev]},{a_time},
                           0,0,nComp,geom[lev],bndry_func,0);
   }
   for (lev = 1; lev <= finest_level; ++lev) {
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy> > crse_bndry_func(geom[lev-1], {m_bcrec_force},
                                                                             PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy> > fine_bndry_func(geom[lev], {m_bcrec_force},
                                                                             PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      Interpolater* mapper = &pc_interp;
      FillPatchTwoLevels(*a_force[lev],IntVect(nGrowForce),a_time,
                         {a_force[lev-1]},{a_time},
                         {a_force[lev]},{a_time},
                         0,0,nComp,geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_force}, 0);
   }
}

#ifdef PLM_USE_EFIELD
void PeleLM::fillpatch_nE(int lev,
                          const amrex::Real a_time,
                          amrex::MultiFab &a_nE,
                          int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirnE>> bndry_func(geom[lev], fetchBCRecArray(NE,1),
                                                                    PeleLMCCFillExtDirnE{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_nE, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->nE),&(m_leveldata_new[lev]->nE)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func, 0);

   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirnE>> crse_bndry_func(geom[lev-1], fetchBCRecArray(NE,1), 
                                                                         PeleLMCCFillExtDirnE{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirnE>> fine_bndry_func(geom[lev], fetchBCRecArray(NE,1),
                                                                         PeleLMCCFillExtDirnE{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_nE, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->nE),&(m_leveldata_new[lev-1]->nE)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->nE),&(m_leveldata_new[lev]->nE)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(NE,1), 0);
   }
}
#endif

#ifdef PLM_USE_EFIELD
void PeleLM::fillpatch_phiV(int lev,
                            const amrex::Real a_time,
                            amrex::MultiFab &a_phiV,
                            int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirPhiV>> bndry_func(geom[lev], fetchBCRecArray(PHIV,1),
                                                                      PeleLMCCFillExtDirPhiV{lprobparm, pmf_data_g, m_nAux});
      FillPatchSingleLevel(a_phiV, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->phiV),&(m_leveldata_new[lev]->phiV)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func, 0);

   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirPhiV>> crse_bndry_func(geom[lev-1], fetchBCRecArray(PHIV,1), 
                                                                           PeleLMCCFillExtDirPhiV{lprobparm, pmf_data_g, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirPhiV>> fine_bndry_func(geom[lev], fetchBCRecArray(PHIV,1),
                                                                           PeleLMCCFillExtDirPhiV{lprobparm, pmf_data_g, m_nAux});
      FillPatchTwoLevels(a_phiV, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->phiV),&(m_leveldata_new[lev-1]->phiV)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->phiV),&(m_leveldata_new[lev]->phiV)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, 1, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(PHIV,1), 0);
   }
}
#endif

// Fill the gradp
void PeleLM::fillpatch_gradp(int lev,
                             const amrex::Real a_time,
                             amrex::MultiFab &a_gp,
                             int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> bndry_func(geom[lev], {m_bcrec_force},
                                                                       PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchSingleLevel(a_gp, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->gp),&(m_leveldata_new[lev]->gp)},
                           {m_t_old[lev], m_t_new[lev]},0,0,AMREX_SPACEDIM,geom[lev], bndry_func, 0);
   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func(geom[lev-1], {m_bcrec_force},
                                                                            PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func(geom[lev], {m_bcrec_force},
                                                                            PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchTwoLevels(a_gp, IntVect(nGhost), a_time,
                         {&(m_leveldata_old[lev-1]->gp),&(m_leveldata_new[lev-1]->gp)},
                         {m_t_old[lev-1], m_t_new[lev-1]},
                         {&(m_leveldata_old[lev]->gp),&(m_leveldata_new[lev]->gp)},
                         {m_t_old[lev], m_t_new[lev]},
                         0, 0, AMREX_SPACEDIM, geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_force}, 0);
   }
}

// Fill the reaction data
void PeleLM::fillpatch_reaction(int lev,
                                const amrex::Real a_time,
                                amrex::MultiFab &a_I_R,
                                int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();
   if (lev == 0) {
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> bndry_func(geom[lev], {m_bcrec_force},
                                                                     PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchSingleLevel(a_I_R, IntVect(nGhost), a_time,
                           {&(m_leveldatareact[lev]->I_R)},{a_time},
                           0, 0, nCompIR(), geom[lev], bndry_func, 0);
   } else {

      // Interpolator
#ifdef AMREX_USE_EB
      Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                             (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
      Interpolater* mapper = &cell_cons_interp;
#endif

      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func(geom[lev-1], {m_bcrec_force},
                                                                            PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func(geom[lev], {m_bcrec_force},
                                                                            PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
      FillPatchTwoLevels(a_I_R, IntVect(nGhost), a_time,
                         {&(m_leveldatareact[lev-1]->I_R)},{a_time},
                         {&(m_leveldatareact[lev]->I_R)},{a_time},
                         0, 0, nCompIR(), geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_force}, 0);
   }
}

// Fill the velocity
void PeleLM::fillcoarsepatch_velocity(int lev,
                                      const amrex::Real a_time,
                                      amrex::MultiFab &a_vel,
                                      int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> crse_bndry_func(geom[lev-1], fetchBCRecArray(VELX,AMREX_SPACEDIM),
                                                                       PeleLMCCFillExtDirVel{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> fine_bndry_func(geom[lev], fetchBCRecArray(VELX,AMREX_SPACEDIM),
                                                                       PeleLMCCFillExtDirVel{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_vel, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->velocity, 0, 0, AMREX_SPACEDIM,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(VELX,AMREX_SPACEDIM), 0);
}

// Fill the mass
void PeleLM::fillcoarsepatch_mass(int lev,
                                  const amrex::Real a_time,
                                  amrex::MultiFab &a_density,
                                  amrex::MultiFab &a_species,
                                  int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   // Density
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDens>> crse_bndry_func_rho(geom[lev-1], fetchBCRecArray(DENSITY,1),
                                                                            PeleLMCCFillExtDirDens{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDens>> fine_bndry_func_rho(geom[lev], fetchBCRecArray(DENSITY,1),
                                                                            PeleLMCCFillExtDirDens{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_density, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->density, 0, 0, 1,
                         geom[lev-1], geom[lev],
                         crse_bndry_func_rho,0,fine_bndry_func_rho,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(DENSITY,1), 0);

   // Species
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirSpec>> crse_bndry_func(geom[lev-1], fetchBCRecArray(FIRSTSPEC,NUM_SPECIES),
                                                                        PeleLMCCFillExtDirSpec{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirSpec>> fine_bndry_func(geom[lev], fetchBCRecArray(FIRSTSPEC,NUM_SPECIES),
                                                                        PeleLMCCFillExtDirSpec{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_species, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->species, 0, 0, NUM_SPECIES,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(FIRSTSPEC,NUM_SPECIES), 0);

}

// Fill the mass
void PeleLM::fillcoarsepatch_energy(int lev,
                                    const amrex::Real a_time,
                                    amrex::MultiFab &a_rhoh,
                                    amrex::MultiFab &a_temp,
                                    int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   // rhoH
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirRhoH>> crse_bndry_func_rhoh(geom[lev-1], fetchBCRecArray(RHOH,1),
                                                                             PeleLMCCFillExtDirRhoH{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirRhoH>> fine_bndry_func_rhoh(geom[lev], fetchBCRecArray(RHOH,1),
                                                                             PeleLMCCFillExtDirRhoH{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_rhoh, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->rhoh, 0, 0, 1,
                         geom[lev-1], geom[lev],
                         crse_bndry_func_rhoh,0,fine_bndry_func_rhoh,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(RHOH,1), 0);

   // temperature
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirTemp>> crse_bndry_func(geom[lev-1], fetchBCRecArray(TEMP,1),
                                                                        PeleLMCCFillExtDirTemp{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirTemp>> fine_bndry_func(geom[lev], fetchBCRecArray(TEMP,1),
                                                                        PeleLMCCFillExtDirTemp{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_temp, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->temp, 0, 0, 1,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(TEMP,1), 0);

}

// Fill the grad P
void PeleLM::fillcoarsepatch_gradp(int lev,
                                   const amrex::Real a_time,
                                   amrex::MultiFab &a_gp,
                                   int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func(geom[lev-1], {m_bcrec_force},
                                                                       PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func(geom[lev], {m_bcrec_force},
                                                                       PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
   InterpFromCoarseLevel(a_gp, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->gp, 0, 0, AMREX_SPACEDIM,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_force}, 0);
}

// Fill the divu
void PeleLM::fillcoarsepatch_divu(int lev,
                                  const amrex::Real a_time,
                                  amrex::MultiFab &a_divu,
                                  int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func(geom[lev-1], {m_bcrec_divu},
                                                                       PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func(geom[lev], {m_bcrec_divu},
                                                                       PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
   InterpFromCoarseLevel(a_divu, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->divu, 0, 0, 1,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_divu}, 0);
}

// Fill coarse patch of reaction
void PeleLM::fillcoarsepatch_reaction(int lev,
                                      const amrex::Real a_time,
                                      amrex::MultiFab &a_I_R,
                                      int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> crse_bndry_func(geom[lev-1], {m_bcrec_force},
                                                                         PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDummy>> fine_bndry_func(geom[lev], {m_bcrec_force},
                                                                         PeleLMCCFillExtDirDummy{lprobparm, m_nAux});
   InterpFromCoarseLevel(a_I_R, IntVect(nGhost), a_time,
                         m_leveldatareact[lev-1]->I_R, 0, 0, nCompIR(),
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, {m_bcrec_force}, 0);
}

#ifdef PLM_USE_EFIELD
// Fill the nE
void PeleLM::fillcoarsepatch_nE(int lev,
                                const amrex::Real a_time,
                                amrex::MultiFab &a_nE,
                                int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirnE>> crse_bndry_func(geom[lev-1], fetchBCRecArray(NE,1),
                                                                      PeleLMCCFillExtDirnE{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirnE>> fine_bndry_func(geom[lev], fetchBCRecArray(NE,1),
                                                                      PeleLMCCFillExtDirnE{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_nE, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->nE, 0, 0, 1,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(NE,1), 0);
}

// Fill the phiV
void PeleLM::fillcoarsepatch_phiV(int lev,
                                  const amrex::Real a_time,
                                  amrex::MultiFab &a_phiV,
                                  int nGhost) {
   ProbParm const* lprobparm = prob_parm.get();

   // Interpolator
#ifdef AMREX_USE_EB
   Interpolater* mapper = (EBFactory(0).isAllRegular()) ?
                          (Interpolater*)(&cell_cons_interp) : (Interpolater*)(&eb_cell_cons_interp);
#else
   Interpolater* mapper = &cell_cons_interp;
#endif

   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirPhiV>> crse_bndry_func(geom[lev-1], fetchBCRecArray(PHIV,1),
                                                                        PeleLMCCFillExtDirPhiV{lprobparm, pmf_data_g, m_nAux});
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirPhiV>> fine_bndry_func(geom[lev], fetchBCRecArray(PHIV,1),
                                                                        PeleLMCCFillExtDirPhiV{lprobparm, pmf_data_g, m_nAux});
   InterpFromCoarseLevel(a_phiV, IntVect(nGhost), a_time,
                         m_leveldata_new[lev-1]->phiV, 0, 0, 1,
                         geom[lev-1], geom[lev],
                         crse_bndry_func,0,fine_bndry_func,0,
                         refRatio(lev-1), mapper, fetchBCRecArray(PHIV,1), 0);
}
#endif

// Fill the physical boundary of a velocity MF
void PeleLM::setPhysBoundaryVel(MultiFab &a_vel,
                                int lev,
                                TimeStamp a_time) {
   BL_PROFILE_VAR("PeleLM::setPhysBoundaryVel()", setPhysBoundaryVel);
   
   Real time = getTime(lev, a_time);

   // Fills interior and periodic domain boundary ghost cells
   a_vel.FillBoundary(geom[lev].periodicity());

   ProbParm const* lprobparm = prob_parm.get();
   PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> bndry_func(geom[lev], fetchBCRecArray(VELX,AMREX_SPACEDIM),
                                                                  PeleLMCCFillExtDirVel{lprobparm, pmf_data_g, m_nAux});

   bndry_func(a_vel,0,AMREX_SPACEDIM,a_vel.nGrowVect(), time, 0);
}
