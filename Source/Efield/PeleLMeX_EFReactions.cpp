#include <PeleLMeX.H>
#include <PeleLMeX_EF_K.H>

using namespace amrex;

void
PeleLM::computeInstantaneousReactionRateEF(
  int lev, const TimeStamp& a_time, MultiFab* a_I_R)
{
  auto ldata_p = getLevelDataPtr(lev, a_time);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(ldata_p->state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    auto const& rhoY = ldata_p->state.const_array(mfi, FIRSTSPEC);
    auto const& rhoH = ldata_p->state.const_array(mfi, RHOH);
    auto const& nE = ldata_p->state.const_array(mfi, NE);
    auto const& T = ldata_p->state.const_array(mfi, TEMP);
    auto const& rhoYdot = a_I_R->array(mfi);
    auto const& nEdot = a_I_R->array(mfi, NUM_SPECIES);

    amrex::ParallelFor(
      bx, [rhoY, rhoH, nE, T, rhoYdot,
           nEdot] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        reactionRateRhoY_EF(i, j, k, rhoY, rhoH, T, nE, rhoYdot, nEdot);
      });
  }
}
