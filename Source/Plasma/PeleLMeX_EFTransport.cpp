#include <PeleLMeX.H>
#include <PeleLMeX_EF_K.H>

void
PeleLM::calcEFTransport(const TimeStamp a_time)
{
  BL_PROFILE("PeleLMeX::calcEFTransport()");

  for (int lev = 0; lev <= finest_level; ++lev) {
    auto ldata_p = getLevelDataPtr(lev, a_time);
    auto dxinv = Geom(lev).InvCellSizeArray();

    auto const& mobE_ma = ldata_p->mobE_cc.arrays();
    auto const& diffE_ma = ldata_p->diffE_cc.arrays();
    auto const& state_ma = ldata_p->state.const_arrays();
    const auto useTab = m_electronKappaTab;
    const auto fixedKe = m_fixedKappaE;

    amrex::ParallelFor(
      ldata_p->diffE_cc, ldata_p->diffE_cc.nGrowVect(),
      [useTab, fixedKe, dxinv, mobE_ma, diffE_ma,
       state_ma] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        constexpr amrex::Real factor = PP_RU_MKS / (Na * elemCharge);
        amrex::Array4<amrex::Real const> rhoY(state_ma[box_no], FIRSTSPEC);
        amrex::Array4<amrex::Real const> phiV(state_ma[box_no], PHIV);
        amrex::Array4<amrex::Real const> T(state_ma[box_no], TEMP);
        getKappaE(
          i, j, k, useTab, fixedKe, dxinv, rhoY, phiV, T, mobE_ma[box_no]);
        getDiffE(i, j, k, factor, T, mobE_ma[box_no], diffE_ma[box_no]);
      });
    amrex::Gpu::streamSynchronize();
  }
}
