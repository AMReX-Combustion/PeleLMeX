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

void PeleLM::getNLStateScaling(Real &nEScale, Real &phiVScale)
{
   Array<Real,2> r = {0.0,0.0};
   for (int comp = 0; comp < 2; comp++) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         if (lev != finest_level) {
            r[comp] = std::max(r[comp],
                               m_leveldatanlsolve[lev]->nlState.norm0(*m_coveredMask[lev],comp,0,true));
         } else {
            r[comp] = std::max(r[comp],
                               m_leveldatanlsolve[lev]->nlState.norm0(comp,0,true,true));
         }
      }
      ParallelDescriptor::ReduceRealMax(r[comp]);
   }
   nEScale = r[0];
   phiVScale = r[1];
}

void PeleLM::getNLResidScaling(Real &nEScale, Real &phiVScale)
{
   Array<Real,2> r = {0.0,0.0};
   for (int comp = 0; comp < 2; comp++) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         if (lev != finest_level) {
            r[comp] = std::max(r[comp], 
                               m_leveldatanlsolve[lev]->nlResid.norm0(*m_coveredMask[lev],comp,0,true));
         } else {
            r[comp] = std::max(r[comp], 
                               m_leveldatanlsolve[lev]->nlResid.norm0(comp,0,true));
         }
      }
      ParallelDescriptor::ReduceRealMax(r[comp]);
   }
   nEScale = r[0];
   phiVScale = r[1];
}

void PeleLM::scaleNLState(const Real &nEScale, const Real &phiVScale)
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      m_leveldatanlsolve[lev]->nlState.mult(1.0/nE_scale,0,1,1);
      m_leveldatanlsolve[lev]->nlState.mult(1.0/phiV_scale,1,1,1);
   }
}

void PeleLM::scaleNLResid(const Vector<MultiFab*> &a_resid, const Real &nEScale, const Real &phiVScale)
{
   for (int lev = 0; lev <= finest_level; ++lev) {
      a_resid[lev]->mult(1.0/FnE_scale,0,1,1);
      a_resid[lev]->mult(1.0/FphiV_scale,1,1,1);
   }
}
