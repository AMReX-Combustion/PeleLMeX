#include <PeleLM.H>

using namespace amrex;

PeleLM::LevelData::LevelData(amrex::BoxArray const& ba,
                             amrex::DistributionMapping const& dm,
                             amrex::FabFactory<FArrayBox> const& factory,
                             int a_incompressible, int a_has_divu, 
                             int a_nAux, int a_nGrowState, int a_nGrowMAC)
{
   velocity.define(ba, dm, AMREX_SPACEDIM, a_nGrowState, MFInfo(), factory);
   gp.define(      ba, dm, AMREX_SPACEDIM, 0           , MFInfo(), factory);
   press.define(   amrex::convert(ba,IntVect::TheNodeVector()),
                       dm, 1             , 1           , MFInfo(), factory);
   visc_cc.define( ba, dm, 1             , 1           , MFInfo(), factory);
   if (! a_incompressible ) {
      density.define (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      species.define (ba, dm, NUM_SPECIES   , a_nGrowState, MFInfo(), factory);
      rhoh.define    (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      rhoRT.define   (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      temp.define    (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      if (a_has_divu) {   
         divu.define (ba, dm, 1             , 1           , MFInfo(), factory);
      }
      diff_cc.define (ba, dm, NUM_SPECIES+2 , 1           , MFInfo(), factory);
#ifdef PLM_USE_EFIELD
      nE.define      (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      phiV.define    (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      diffE_cc.define(ba, dm, 1             , 1           , MFInfo(), factory);
      mobE_cc.define (ba, dm, 1             , 1           , MFInfo(), factory);
      mob_cc.define  (ba, dm, NUM_SPECIES   , 1           , MFInfo(), factory);
#endif
   }
   if ( a_nAux > 0 ) {
      auxiliaries.define(ba, dm, a_nAux, a_nGrowState, MFInfo(), factory);
   }
}

PeleLM::LevelDataReact::LevelDataReact(const amrex::BoxArray &ba,
                                       const amrex::DistributionMapping &dm, 
                                       const amrex::FabFactory<FArrayBox> &factory)
{
   I_R.define(ba, dm, NUM_SPECIES, 0, MFInfo(), factory);
   functC.define(ba, dm, 1, 0, MFInfo(), factory);
#ifdef PLM_USE_EFIELD
   I_RnE.define(ba, dm, 1, 0, MFInfo(), factory);
#endif
}

#ifdef PLM_USE_EFIELD
PeleLM::LevelDataNLSolve::LevelDataNLSolve(amrex::BoxArray const& ba,
                                           amrex::DistributionMapping const& dm,
                                           amrex::FabFactory<FArrayBox> const& factory,
                                           int a_nGrow)
{
   nlState.define(ba, dm, 2, a_nGrow, MFInfo(), factory);
   nlResid.define(ba, dm, 2, a_nGrow, MFInfo(), factory);
   backgroundCharge.define(ba, dm, 1, 0, MFInfo(), factory);
   for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      const BoxArray& faceba = amrex::convert(ba,IntVect::TheDimensionVector(idim));
      gPhiVOld[idim].define(faceba, dm, 1, 0, MFInfo(), factory);
      uEffnE[idim].define(faceba, dm, 1, 0, MFInfo(), factory);
      umac[idim].define(faceba, dm, 1, 0, MFInfo(), factory);
   }
}
#endif

PeleLM::AdvanceDiffData::AdvanceDiffData(int a_finestLevel,
                                         const amrex::Vector<amrex::BoxArray> &ba,
                                         const amrex::Vector<amrex::DistributionMapping> &dm, 
                                         const amrex::Vector<std::unique_ptr<amrex::FabFactory<FArrayBox>>> &factory,
                                         int nGrowAdv,
                                         int a_use_wbar,
                                         int is_init)
{
   if (is_init) {                   // All I need is a container for a single diffusion term
      // Resize Vectors
      Dnp1.resize(a_finestLevel+1);

      // Define MFs
      for (int lev = 0; lev <= a_finestLevel; lev++ ) {
         Dnp1[lev].define(ba[lev], dm[lev], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[lev]);
      }
   } else {
      // Resize Vectors
      Dn.resize(a_finestLevel+1);
      Dnp1.resize(a_finestLevel+1);
      Dhat.resize(a_finestLevel+1);
      if ( a_use_wbar ) {
         Dwbar.resize(a_finestLevel+1);
         wbar_fluxes.resize(a_finestLevel+1);
      }

      // Define MFs
      for (int lev = 0; lev <= a_finestLevel; lev++ ) {
         Dn[lev].define(ba[lev], dm[lev], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[lev]);
         Dnp1[lev].define(ba[lev], dm[lev], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[lev]);
         Dhat[lev].define(ba[lev], dm[lev], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[lev]);
         if ( a_use_wbar ) {
            Dwbar[lev].define(ba[lev], dm[lev], NUM_SPECIES, nGrowAdv, MFInfo(), *factory[lev]);
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
               const BoxArray& faceba = amrex::convert(ba[lev],IntVect::TheDimensionVector(idim));
               wbar_fluxes[lev][idim].define(faceba,dm[lev], NUM_SPECIES, 0, MFInfo(), *factory[lev]);
            }
         }
      }
   }
}

PeleLM::AdvanceAdvData::AdvanceAdvData(int a_finestLevel,
                                       const amrex::Vector<amrex::BoxArray> &ba,
                                       const amrex::Vector<amrex::DistributionMapping> &dm, 
                                       const amrex::Vector<std::unique_ptr<amrex::FabFactory<FArrayBox>>> &factory,
                                       int a_incompressible,
                                       int nGrowAdv,
                                       int nGrowMAC)
{
   // Resize Vectors
   umac.resize(a_finestLevel+1);
   AofS.resize(a_finestLevel+1);
   if ( !a_incompressible ) {
      chi.resize(a_finestLevel+1);
      Forcing.resize(a_finestLevel+1);
      mac_divu.resize(a_finestLevel+1);
   }

   // Define MFs
   for (int lev = 0; lev <= a_finestLevel; lev++ ) {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
         const BoxArray& faceba = amrex::convert(ba[lev],IntVect::TheDimensionVector(idim));
         umac[lev][idim].define(faceba,dm[lev], 1, nGrowMAC, MFInfo(), *factory[lev]);
      }
      if ( a_incompressible ) {
         AofS[lev].define(ba[lev], dm[lev], AMREX_SPACEDIM , 0, MFInfo(), *factory[lev]);
      } else {
         AofS[lev].define(ba[lev], dm[lev], NVAR , 0, MFInfo(), *factory[lev]);
         chi[lev].define(ba[lev], dm[lev], 1, 1, MFInfo(), *factory[lev]);
         Forcing[lev].define(ba[lev], dm[lev], NUM_SPECIES+1, nGrowAdv, MFInfo(), *factory[lev]); // Species + TEMP
         mac_divu[lev].define(ba[lev], dm[lev], 1, nGrowAdv, MFInfo(), *factory[lev]);
      }
   }
}

void
PeleLM::copyStateNewToOld(int nGhost) {
   AMREX_ASSERT(nGhost<=m_nGrowState);
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_old[lev]->velocity,m_leveldata_new[lev]->velocity,0,0,AMREX_SPACEDIM,nGhost);
      if ( !m_incompressible ) {
         MultiFab::Copy(m_leveldata_old[lev]->density,m_leveldata_new[lev]->density,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->species,m_leveldata_new[lev]->species,0,0,NUM_SPECIES,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->rhoh,m_leveldata_new[lev]->rhoh,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->temp,m_leveldata_new[lev]->temp,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->rhoRT,m_leveldata_new[lev]->rhoRT,0,0,1,nGhost);
         if ( m_has_divu ) {
            MultiFab::Copy(m_leveldata_old[lev]->divu,m_leveldata_new[lev]->divu,0,0,1,std::min(nGhost,1));
         }
#ifdef PLM_USE_EFIELD
         MultiFab::Copy(m_leveldata_old[lev]->nE,m_leveldata_new[lev]->nE,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->phiV,m_leveldata_new[lev]->phiV,0,0,1,nGhost);
#endif
      }
   }
}

void
PeleLM::copyPressNewToOld() {
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_old[lev]->press,m_leveldata_new[lev]->press,0,0,1,1);
      MultiFab::Copy(m_leveldata_old[lev]->gp,m_leveldata_new[lev]->gp,0,0,AMREX_SPACEDIM,0);
   }
}

void
PeleLM::copyStateOldToNew(int nGhost) {
   AMREX_ASSERT(nGhost<=m_nGrowState);
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_new[lev]->velocity,m_leveldata_old[lev]->velocity,0,0,AMREX_SPACEDIM,nGhost);
      if ( !m_incompressible ) {
         MultiFab::Copy(m_leveldata_new[lev]->density,m_leveldata_old[lev]->density,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->species,m_leveldata_old[lev]->species,0,0,NUM_SPECIES,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->rhoh,m_leveldata_old[lev]->rhoh,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->temp,m_leveldata_old[lev]->temp,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->rhoRT,m_leveldata_old[lev]->rhoRT,0,0,1,nGhost);
         if ( m_has_divu ) {
            MultiFab::Copy(m_leveldata_new[lev]->divu,m_leveldata_old[lev]->divu,0,0,1,std::min(nGhost,1));
         }
#ifdef PLM_USE_EFIELD
         MultiFab::Copy(m_leveldata_new[lev]->nE,m_leveldata_old[lev]->nE,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->phiV,m_leveldata_old[lev]->phiV,0,0,1,nGhost);
#endif
      }
   }
}


void
PeleLM::copyTransportOldToNew() {
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_new[lev]->visc_cc,m_leveldata_old[lev]->visc_cc,0,0,1,1);
      if ( !m_incompressible ) {
         MultiFab::Copy(m_leveldata_new[lev]->diff_cc,m_leveldata_old[lev]->diff_cc,0,0,NUM_SPECIES+2,1);
#ifdef PLM_USE_EFIELD
         MultiFab::Copy(m_leveldata_new[lev]->diffE_cc,m_leveldata_old[lev]->diffE_cc,0,0,1,1);
         MultiFab::Copy(m_leveldata_new[lev]->mobE_cc,m_leveldata_old[lev]->mobE_cc,0,0,1,1);
         MultiFab::Copy(m_leveldata_new[lev]->mob_cc,m_leveldata_old[lev]->mob_cc,0,0,NUM_SPECIES,1);
#endif
      }
   }
}

void
PeleLM::copyDiffusionOldToNew(std::unique_ptr<AdvanceDiffData> &diffData) {
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(diffData->Dnp1[lev],diffData->Dn[lev],0,0,NUM_SPECIES+2,m_nGrowAdv);
   }
}
