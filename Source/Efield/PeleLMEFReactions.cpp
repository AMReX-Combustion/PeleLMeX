#include <PeleLM.H>
#include <PeleLMEF_K.H>

using namespace amrex;

void PeleLM::computeInstantaneousReactionRateEF(int lev,
                                                const TimeStamp &a_time,
                                                MultiFab* a_I_R)
{
   auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldata_p->species,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhoY    = ldata_p->species.const_array(mfi);
      auto const& rhoH    = ldata_p->rhoh.const_array(mfi);
      auto const& nE      = ldata_p->nE.const_array(mfi);
      auto const& T       = ldata_p->temp.const_array(mfi);
      auto const& rhoYdot = a_I_R->array(mfi);
      auto const& nEdot   = a_I_R->array(mfi,NUM_SPECIES);

      amrex::ParallelFor(bx, [rhoY, rhoH, nE, T, rhoYdot, nEdot]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         reactionRateRhoY_EF( i, j, k, rhoY, rhoH, T, nE, rhoYdot, nEdot );
      });
   }
}
