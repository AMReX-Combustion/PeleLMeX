#include <PeleLMeX.H>
#include <PeleLMeX_EF_K.H>

void
PeleLM::calcEFTransport(const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::calcEFTransport()");

  for (int lev = 0; lev <= finest_level; ++lev) {
    auto ldata_p = getLevelDataPtr(lev, a_time);
    auto dxinv = Geom(lev).InvCellSizeArray();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(ldata_p->diffE_cc, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      const amrex::Box& gbx = mfi.growntilebox();
      auto const& mobE = ldata_p->mobE_cc.array(mfi);
      auto const& diffE = ldata_p->diffE_cc.array(mfi);
      auto const& rhoY = ldata_p->state.const_array(mfi, FIRSTSPEC);
      auto const& phiV = ldata_p->state.const_array(mfi, PHIV);
      auto const& T = ldata_p->state.const_array(mfi, TEMP);
      amrex::Real factor = PP_RU_MKS / (Na * elemCharge);
      const auto useTab = m_electronKappaTab;
      const auto fixedKe = m_fixedKappaE;
      amrex::ParallelFor(
        gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          getKappaE(i, j, k, useTab, fixedKe, dxinv, rhoY, phiV, T, mobE);
          getDiffE(i, j, k, factor, T, mobE, diffE);
        });
    }
  }
}
