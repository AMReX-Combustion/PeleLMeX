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
      density.define(ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      species.define(ba, dm, NUM_SPECIES   , a_nGrowState, MFInfo(), factory);
      rhoh.define   (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      rhoRT.define  (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      temp.define   (ba, dm, 1             , a_nGrowState, MFInfo(), factory);
      if (a_has_divu) {   
         divu.define(ba, dm, 1             , 1           , MFInfo(), factory);
      }
      diff_cc.define(ba, dm, NUM_SPECIES+2 , 1           , MFInfo(), factory);
   }
   if ( a_nAux > 0 ) {
      auxiliaries.define(ba, dm, a_nAux, a_nGrowState, MFInfo(), factory);
   }
}

PeleLM::AdvanceDiffData::AdvanceDiffData(int a_finestLevel,
                                         const amrex::Vector<amrex::BoxArray> &ba,
                                         const amrex::Vector<amrex::DistributionMapping> &dm, 
                                         const amrex::Vector<std::unique_ptr<amrex::FabFactory<FArrayBox>>> &factory,
                                         int nGrowAdv,
                                         int a_use_wbar)
{
   // Resize Vectors
   Dn.resize(a_finestLevel+1);
   Dnp1.resize(a_finestLevel+1);
   Dhat.resize(a_finestLevel+1);
   if ( a_use_wbar ) {
      Dwbar.resize(a_finestLevel+1);
   }

   // Define MFs
   for (int n = 0; n <= a_finestLevel; n++ ) {
      Dn[n].define(ba[n], dm[n], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[n]);
      Dnp1[n].define(ba[n], dm[n], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[n]);
      Dhat[n].define(ba[n], dm[n], NUM_SPECIES+2 , nGrowAdv, MFInfo(), *factory[n]);
      if ( a_use_wbar ) {
         Dwbar[n].define(ba[n], dm[n], NUM_SPECIES, nGrowAdv, MFInfo(), *factory[n]);
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
         Forcing[lev].define(ba[lev], dm[lev], NUM_SPECIES+1, nGrowAdv, MFInfo(), *factory[lev]); // Species + RHOH
      }
   }
}

void
PeleLM::copyStateNewToOld(int nGhost) {
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_old[lev]->velocity,m_leveldata_new[lev]->velocity,0,0,AMREX_SPACEDIM,nGhost);
      if ( !m_incompressible ) {
         MultiFab::Copy(m_leveldata_old[lev]->density,m_leveldata_new[lev]->density,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->species,m_leveldata_new[lev]->species,0,0,NUM_SPECIES,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->rhoh,m_leveldata_new[lev]->rhoh,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->temp,m_leveldata_new[lev]->temp,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_old[lev]->rhoRT,m_leveldata_new[lev]->rhoRT,0,0,1,nGhost);
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
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_new[lev]->velocity,m_leveldata_old[lev]->velocity,0,0,AMREX_SPACEDIM,nGhost);
      if ( !m_incompressible ) {
         MultiFab::Copy(m_leveldata_new[lev]->density,m_leveldata_old[lev]->density,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->species,m_leveldata_old[lev]->species,0,0,NUM_SPECIES,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->rhoh,m_leveldata_old[lev]->rhoh,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->temp,m_leveldata_old[lev]->temp,0,0,1,nGhost);
         MultiFab::Copy(m_leveldata_new[lev]->rhoRT,m_leveldata_old[lev]->rhoRT,0,0,1,nGhost);
      }
   }
}


void
PeleLM::copyTransportOldToNew() {
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(m_leveldata_new[lev]->visc_cc,m_leveldata_old[lev]->visc_cc,0,0,1,1);
      if ( !m_incompressible ) {
         MultiFab::Copy(m_leveldata_new[lev]->diff_cc,m_leveldata_old[lev]->diff_cc,0,0,NUM_SPECIES+2,1);
      }
   }
}

void
PeleLM::copyDiffusionOldToNew(std::unique_ptr<AdvanceDiffData> &diffData) {
   for (int lev = 0; lev <= finest_level; lev++ ) {
      MultiFab::Copy(diffData->Dnp1[lev],diffData->Dn[lev],0,0,NUM_SPECIES+2,m_nGrowAdv);
   }
}
