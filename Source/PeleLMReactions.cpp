#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

void PeleLM::advanceChemistry()
{
}

void PeleLM::computeInstantaneousReactionRate(const Vector<MultiFab*> &I_R,
                                              const TimeStamp &a_time)
{
   BL_PROFILE_VAR("PeleLM::computeInstantaneousReactionRate()", computeInstantaneousReactionRate);

   for (int lev = 0; lev <= finest_level; ++lev) {
      // Setup covered mask
      MultiFab mask(grids[lev],dmap[lev],1,0);
      mask.setVal(1.0);
      computeInstantaneousReactionRate(lev, a_time, mask, I_R[lev]);
   }
}

void PeleLM::computeInstantaneousReactionRate(int lev,
                                              const TimeStamp &a_time,
                                              const MultiFab &a_mask,
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
      auto const& T       = ldata_p->temp.const_array(mfi);
      auto const& mask    = a_mask.const_array(mfi);
      auto const& rhoYdot = a_I_R->array(mfi);

      amrex::ParallelFor(bx, [rhoY, rhoH, T, mask, rhoYdot]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {    
         reactionRateRhoY( i, j, k, rhoY, rhoH, T, mask,
                           rhoYdot );
      });  
   }
}
