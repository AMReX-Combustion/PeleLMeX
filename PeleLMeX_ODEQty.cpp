#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

#if NUM_ODE > 0
void
PeleLM::predictODEQty()
{
  // Uses forward Euler to predict values for ODE qty at tnp1
  // If m_ext_sources_SDC = false, no SDC corrections used
  const auto dt = m_dt;
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto const& state_arrs = getLevelDataPtr(lev, AmrNewTime)->state.arrays();
    auto const& ext_src_arrs = m_extSource[lev]->arrays();
    amrex::ParallelFor(
      *m_extSource[lev], amrex::IntVect(0), NUM_ODE,
      [state_arrs, ext_src_arrs,
       dt] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
        const amrex::Real S_ext_n = ext_src_arrs[box_no](i, j, k, FIRSTODE + n);
        state_arrs[box_no](i, j, k, FIRSTODE + n) += dt * S_ext_n;
      });
    amrex::Gpu::streamSynchronize();
  }
}
#endif
