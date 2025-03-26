#include <PeleLMeX.H>
#include <PeleLMeX_EF_K.H>

using namespace amrex;

void
PeleLM::calcEFTransport(const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::calcEFTransport()");

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto ldata_p = getLevelDataPtr(lev, a_time);
    auto dxinv = Geom(lev).InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldata_p->diffE_cc, TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const Box& gbx = mfi.growntilebox();
      auto const& mobE = ldata_p->mobE_cc.array(mfi);
      auto const& diffE = ldata_p->diffE_cc.array(mfi);
      auto const& rhoY = ldata_p->state.const_array(mfi, FIRSTSPEC);
      auto const& phiV = ldata_p->state.const_array(mfi, PHIV);
      auto const& T = ldata_p->state.const_array(mfi, TEMP);
      Real factor = PP_RU_MKS / (Na * elemCharge);
      amrex::ParallelFor(
        gbx,
        [mobE, diffE, rhoY, phiV, T, factor, dxinv, useTab = m_electronKappaTab,
         fixedKe =
           m_fixedKappaE] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          getKappaE(i, j, k, useTab, fixedKe, dxinv, rhoY, phiV, T, mobE);
          getDiffE(i, j, k, factor, T, mobE, diffE);
        });
    }
  }
}
