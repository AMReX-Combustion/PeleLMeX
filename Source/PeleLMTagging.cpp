#include <PeleLM.H>
#include <AMReX_TagBox.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

using namespace amrex;

void
PeleLM::ErrorEst( int lev,
                  TagBoxArray& tags,
                  Real time,
                  int ng)
{
   BL_PROFILE_VAR("PeleLM::ErrorEst()", ErrorEst);

#ifdef AMREX_USE_EB
   if (m_refine_cutcells) {
      const MultiFab& rho = (getLevelDataPtr(lev,AmrNewTime))->density;
      TagCutCells(tags, rho);
   }
#endif

   for (int n = 0; n < errTags.size(); ++n) {
      std::unique_ptr<MultiFab> mf;
      if (errTags[n].Field() != std::string()) {
         mf = deriveComp(errTags[n].Field(), time, lev, errTags[n].NGrow());
      }
      errTags[n](tags,mf.get(),TagBox::CLEAR,TagBox::SET,time,lev,geom[lev]);
   }
}
