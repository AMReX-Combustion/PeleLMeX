#include <PeleLM.H>
#include <PeleLMEF_K.H>

using namespace amrex;

void PeleLM::calcEFTransport(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::calcEFTransport()", calcEFTransport);

   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->diffE_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx     = mfi.growntilebox();
         auto const& mobE   = ldata_p->mobE_cc.array(mfi);
         auto const& diffE  = ldata_p->diffE_cc.array(mfi);
         auto const& T      = ldata_p->state.const_array(mfi,TEMP);
         Real factor = PP_RU_MKS / ( Na * elemCharge );
         amrex::ParallelFor(gbx, [mobE, diffE, T, factor]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getKappaE(i,j,k,mobE);
            getDiffE(i,j,k,factor,T,mobE,diffE);
         });
      }
   }
}
