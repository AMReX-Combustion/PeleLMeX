#include <PeleLM.H>
#include <PeleLM_K.H>
#include <PeleLMEF_K.H>

using namespace amrex;

Vector<Array<MultiFab*,AMREX_SPACEDIM>>
PeleLM::getNLgradPhiVVect() {
   Vector<Array<MultiFab*,AMREX_SPACEDIM>> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(GetArrOfPtrs(m_leveldatanlsolve[lev]->gPhiVOld));
   }
   return r;
}

Vector<Array<MultiFab*,AMREX_SPACEDIM>>
PeleLM::getUeffVect() {
   Vector<Array<MultiFab*,AMREX_SPACEDIM>> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(GetArrOfPtrs(m_leveldatanlsolve[lev]->uEffnE));
   }
   return r;
}

Vector<MultiFab*>
PeleLM::getNLresidVect() {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldatanlsolve[lev]->nlResid));
   }
   return r;
}

Vector<MultiFab*>
PeleLM::getNLstateVect() {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      r.push_back(&(m_leveldatanlsolve[lev]->nlState));
   }
   return r;
}
