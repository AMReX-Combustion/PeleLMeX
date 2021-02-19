#include <PeleLM.H>

using namespace amrex;

PeleLM::PeleLM() = default;

PeleLM::~PeleLM() = default;

Vector<MultiFab *>
PeleLM::getSpeciesVect(TimeStamp a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->species));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->species));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDensityVect(TimeStamp a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->density));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->density));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getTempVect(TimeStamp a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->temp));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->temp));
      }
   }
   return r;
}

Vector<MultiFab *>
PeleLM::getDiffusivityVect(TimeStamp a_time) {
   Vector<MultiFab*> r;
   r.reserve(finest_level+1);
   if ( a_time == AmrOldTime ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_old[lev]->diff_cc));
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         r.push_back(&(m_leveldata_new[lev]->diff_cc));
      }
   }
   return r;
}
