#include <PeleLMeX.H>
#include <PeleLMeX_EF_K.H>

void
PeleLM::computeInstantaneousReactionRateEF(
  const int lev, const TimeStamp a_time, amrex::MultiFab* a_I_R)
{
  auto ldata_p = getLevelDataPtr(lev, a_time);

  auto const& state_ma = ldata_p->state.const_arrays();
  auto const& I_R_ma = a_I_R->arrays();

  amrex::ParallelFor(
    ldata_p->state, [state_ma, I_R_ma] AMREX_GPU_DEVICE(
                      int box_no, int i, int j, int k) noexcept {
      amrex::Array4<amrex::Real const> rhoY(state_ma[box_no], FIRSTSPEC);
      amrex::Array4<amrex::Real const> rhoH(state_ma[box_no], RHOH);
      amrex::Array4<amrex::Real const> nE(state_ma[box_no], NE);
      amrex::Array4<amrex::Real const> T(state_ma[box_no], TEMP);
      amrex::Array4<amrex::Real> rhoYdot(I_R_ma[box_no], 0);
      amrex::Array4<amrex::Real> nEdot(I_R_ma[box_no], NUM_SPECIES);
      reactionRateRhoY_EF(i, j, k, rhoY, rhoH, T, nE, rhoYdot, nEdot);
    });
  amrex::Gpu::streamSynchronize();
}
