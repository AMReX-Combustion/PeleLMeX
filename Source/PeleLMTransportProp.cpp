#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

void PeleLM::calcDiffusivity(TimeStamp time) {
   BL_PROFILE_VAR("PeleLM::calcDiffusivity()", calcDiffusivity);

   for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldata_p = getLevelDataPtr(lev,time); 

      // Transport data pointer
      TransParm const* ltransparm = trans_parm_g;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx     = mfi.growntilebox();
         auto const& rhoY   = ldata_p->species.const_array(mfi);
         auto const& T      = ldata_p->temp.const_array(mfi);
         auto const& rhoD   = ldata_p->diff_cc.array(mfi,0);
         auto const& lambda = ldata_p->diff_cc.array(mfi,NUM_SPECIES);
         auto const& mu     = ldata_p->diff_cc.array(mfi,NUM_SPECIES+1);

         // TODO: unity Lewis

         amrex::ParallelFor(gbx, [rhoY, T, rhoD, lambda, mu, ltransparm]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getTransportCoeff( i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
         });
      }
   }
}
