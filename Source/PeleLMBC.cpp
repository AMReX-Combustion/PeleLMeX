#include <PeleLM.H>
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


void PeleLM::setBoundaryConditions() {

   // Initialize the state BCRec
   m_bcrec_state.resize(NVAR);

   // Convert m_phys_bc into field BCs
   // Get m_phys_bc
   const int* lo_bc = m_phys_bc.lo();
   const int* hi_bc = m_phys_bc.hi();

   // Velocity
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      for (int idim2 = 0; idim2 < AMREX_SPACEDIM; idim2++) {
         if ( idim == idim2 ) {
            m_bcrec_state[VELX+idim].setLo(idim2,norm_vel_bc[lo_bc[idim2]]);
            m_bcrec_state[VELX+idim].setHi(idim2,norm_vel_bc[lo_bc[idim2]]);
         } else {
            m_bcrec_state[VELX+idim].setLo(idim2,tang_vel_bc[lo_bc[idim2]]);
            m_bcrec_state[VELX+idim].setHi(idim2,tang_vel_bc[lo_bc[idim2]]);
         }
      }
   }  

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

// Fill the entire state at once
void PeleLM::fillPatchState(int lev, TimeStamp a_time) {
   BL_PROFILE_VAR("PeleLM::fillPatchState()", fillPatchState);
   
   auto ldata_p = getLevelDataPtr(lev,a_time);
   Real time = getTime(lev, a_time);

   fillpatch_velocity(lev, time, ldata_p->velocity, m_nGrowState);
   fillpatch_mass(lev, time, ldata_p->density, ldata_p->species, m_nGrowState);
   fillpatch_energy(lev, time, ldata_p->rhoh, ldata_p->temp, m_nGrowState);
}

// Fill the velocity
void PeleLM::fillpatch_velocity(int lev,
                                const amrex::Real a_time,
                                amrex::MultiFab &a_vel,
                                int nGhost) { 
   if (lev == 0) {
      ProbParm const* lprobparm = prob_parm.get();
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirVel>> bndry_func(geom[lev], fetchBCRecArray(VELX,AMREX_SPACEDIM), 
                                                                     PeleLMCCFillExtDirVel{lprobparm, m_nAux});
      FillPatchSingleLevel(a_vel, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->velocity),&(m_leveldata_new[lev]->velocity)},
                           {m_t_old[lev], m_t_new[lev]},0,0,AMREX_SPACEDIM,geom[lev], bndry_func, 0);
   } else {
      Abort("FillPatchtwoLevels for velocity need to be implemented");
   }
}

// Fill the density and mass fractions only
void PeleLM::fillpatch_mass(int lev,
                            const amrex::Real a_time,
                            amrex::MultiFab &a_density,
                            amrex::MultiFab &a_species,
                            int nGhost) { 
   if (lev == 0) {
      ProbParm const* lprobparm = prob_parm.get();

      // Density
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirDens>> bndry_func_rho(geom[lev], fetchBCRecArray(DENSITY,1), 
                                                                          PeleLMCCFillExtDirDens{lprobparm, m_nAux});
      FillPatchSingleLevel(a_density, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->density),&(m_leveldata_new[lev]->density)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func_rho, 0);

      // Species
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirSpec>> bndry_func(geom[lev], fetchBCRecArray(FIRSTSPEC,NUM_SPECIES), 
                                                                      PeleLMCCFillExtDirSpec{lprobparm, m_nAux});
      FillPatchSingleLevel(a_species, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->species),&(m_leveldata_new[lev]->species)},
                           {m_t_old[lev], m_t_new[lev]},0,0,NUM_SPECIES,geom[lev], bndry_func, 0);
   } else {
      Abort("FillPatchtwoLevels for mass need to be implemented");
   }
}

// Fill rhoH and temperature
void PeleLM::fillpatch_energy(int lev,
                              const amrex::Real a_time,
                              amrex::MultiFab &a_rhoh,
                              amrex::MultiFab &a_temp,
                              int nGhost) { 
   if (lev == 0) {
      ProbParm const* lprobparm = prob_parm.get();

      // rhoH
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirRhoH>> bndry_func_rhoh(geom[lev], fetchBCRecArray(RHOH,1), 
                                                                          PeleLMCCFillExtDirRhoH{lprobparm, m_nAux});
      FillPatchSingleLevel(a_rhoh, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->rhoh),&(m_leveldata_new[lev]->rhoh)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func_rhoh, 0);

      // temperature
      PhysBCFunct<GpuBndryFuncFab<PeleLMCCFillExtDirTemp>> bndry_func(geom[lev], fetchBCRecArray(TEMP,1), 
                                                                      PeleLMCCFillExtDirTemp{lprobparm, m_nAux});
      FillPatchSingleLevel(a_temp, IntVect(nGhost), a_time,
                           {&(m_leveldata_old[lev]->temp),&(m_leveldata_new[lev]->temp)},
                           {m_t_old[lev], m_t_new[lev]},0,0,1,geom[lev], bndry_func, 0);
   } else {
      Abort("FillPatchtwoLevels for energy need to be implemented");
   }
}

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
                                                                  PeleLMCCFillExtDirVel{lprobparm, m_nAux});

   bndry_func(a_vel,0,AMREX_SPACEDIM,a_vel.nGrowVect(), time, 0);
}

