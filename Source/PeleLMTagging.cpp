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
                  int /*ng*/)
{
   BL_PROFILE_VAR("PeleLM::ErrorEst()", ErrorEst);

#ifdef AMREX_USE_EB
   // Tag EB up to m_EB_refine_LevMax-1 if Static or
   //              m_EB_refine_LevAdapt-1 if Adaptive
   if ( ( m_EB_refine_type == "Static" && lev < m_EB_refine_LevMax ) ||
        ( m_EB_refine_type == "Adaptive" && lev < m_EB_refine_LevAdapt ) ) {
      const MultiFab& state = (getLevelDataPtr(lev,AmrNewTime))->state;
      TagCutCells(tags, state);
   }
#endif

   for (int n = 0; n < errTags.size(); ++n) {
      std::unique_ptr<MultiFab> mf;
      if (errTags[n].Field() != std::string()) {
         mf = deriveComp(errTags[n].Field(), time, lev, errTags[n].NGrow());
      }
      errTags[n](tags,mf.get(),TagBox::CLEAR,TagBox::SET,time,lev,geom[lev]);
   }

#ifdef AMREX_USE_EB
   // Untag cell close to EB
   if ( m_EB_refine_type == "Static" && lev >= m_EB_refine_LevMax ) {
      // Get distance function at current level
      MultiFab signDist(grids[lev],dmap[lev],1,0,MFInfo(),EBFactory(lev));
      getEBDistance(lev, signDist);
      //VisMF::Write(signDist,"signDistLev"+std::to_string(lev));

      // Estimate how far I need to derefine
      Real diagFac = std::sqrt(2.0) * m_derefineEBBuffer;
      Real clearTagDist = Geom(m_EB_refine_LevMax).CellSize(0) * static_cast<Real>(nErrorBuf(m_EB_refine_LevMax)) * diagFac;
      for (int ilev = m_EB_refine_LevMax+1; ilev <= finest_level; ++ilev) {
          clearTagDist += static_cast<Real>(nErrorBuf(ilev)) * Geom(m_EB_refine_LevMax).CellSize(0) * diagFac;
      }
      //Print() << " clearTagDist " <<  clearTagDist << "\n";

      // Untag cells too close to EB
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(tags,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const auto& bx    = mfi.tilebox();
          const auto& dist  = signDist.const_array(mfi);
          auto tag          = tags.array(mfi);
          amrex::ParallelFor(bx,
          [=] AMREX_GPU_HOST_DEVICE (int i, int j, int k)
          {
              if (dist(i,j,k) < clearTagDist ) {
                  tag(i,j,k) = TagBox::CLEAR;
              }
          });
      }
   }
#endif
}
