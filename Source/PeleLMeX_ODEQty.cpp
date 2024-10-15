#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

using namespace amrex;

void
PeleLM::predictODEQty()
{
  // Uses forward Euler to predict values for ODE qty at tnp1
  // If m_ext_sources_SDC = false, no SDC corrections used
  for (int lev = 0; lev <= finest_level; lev++) {
    auto const& state_arrs = getLevelDataPtr(lev, AmrNewTime)->state.arrays();
    auto const& ext_src_arrs = m_extSource[lev]->arrays();
    ParallelFor(
      *m_extSource[lev], [state_arrs, ext_src_arrs, dt = m_dt] AMREX_GPU_DEVICE(
                           int box_no, int i, int j, int k) noexcept {
        for (int n = 0; n < NUM_ODE; n++) {
          Real const& B_n = state_arrs[box_no](i, j, k, FIRSTODE + n);
          Real const& S_ext_n = ext_src_arrs[box_no](i, j, k, FIRSTODE + n);
          state_arrs[box_no](i, j, k, FIRSTODE + n) = B_n + dt * S_ext_n;
        }
      });
    Gpu::streamSynchronize();
  }
}