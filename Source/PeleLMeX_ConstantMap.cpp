#include <PeleLMeX_ConstantMap.H>

#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>
#include <limits>

using amrex::Real;

ConstantMap::ConstantMap()
{
  amrex::ParmParse pp("ConstantMap");
  m_fac.assign(AMREX_SPACEDIM, 1.0);
  pp.queryarr("scaling_factor", m_fac, 0, AMREX_SPACEDIM);

  // Validate that all scaling factors are strictly positive
  // to avoid division-by-zero or NaN issues in downstream operations
  // (MAC projection, diffusion, etc. divide by fac and fac^2)
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (m_fac[d] <= std::numeric_limits<amrex::Real>::epsilon()) {
      amrex::Abort(
        "ConstantMap: All scaling factors must be strictly positive.");
    }
  }
}

void
ConstantMap::create_map(int lev, const amrex::Geometry& /*geom*/)
{
  AMREX_ASSERT(lev >= 0 && lev < num_levels());
  fill_cc(lev);
  fill_nd(lev);
  fill_fc(lev);
  amrex::Gpu::streamSynchronize();
}

MeshMapEvaluator
ConstantMap::make_evaluator() const
{
  MeshMapEvaluator e;
  e.m_kind = MeshMapEvaluator::Kind::Constant;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    e.m_p[d] = m_fac[d];
    e.m_q[d] = -1;
  }
  return e;
}

void
ConstantMap::fill_cc(int lev)
{
  AMREX_D_TERM(
    const Real fx = m_fac[0];, const Real fy = m_fac[1];
    , const Real fz = m_fac[2];);
  const Real detJ = AMREX_D_TERM(fx, *fy, *fz);

  const auto& fac_ma = m_fac_cc[lev].arrays();
  const auto& det_ma = m_detJ_cc[lev].arrays();
  amrex::ParallelFor(
    m_fac_cc[lev], m_fac_cc[lev].nGrowVect(),
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      auto f = fac_ma[box_no];
      auto d = det_ma[box_no];
      AMREX_D_TERM(f(i, j, k, 0) = fx;, f(i, j, k, 1) = fy;
                   , f(i, j, k, 2) = fz;);
      d(i, j, k) = detJ;
    });
}

void
ConstantMap::fill_nd(int lev)
{
  AMREX_D_TERM(
    const Real fx = m_fac[0];, const Real fy = m_fac[1];
    , const Real fz = m_fac[2];);
  const Real detJ = AMREX_D_TERM(fx, *fy, *fz);

  const auto& fac_ma = m_fac_nd[lev].arrays();
  const auto& det_ma = m_detJ_nd[lev].arrays();
  amrex::ParallelFor(
    m_fac_nd[lev], m_fac_nd[lev].nGrowVect(),
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      auto f = fac_ma[box_no];
      auto d = det_ma[box_no];
      AMREX_D_TERM(f(i, j, k, 0) = fx;, f(i, j, k, 1) = fy;
                   , f(i, j, k, 2) = fz;);
      d(i, j, k) = detJ;
    });
}

void
ConstantMap::fill_fc(int lev)
{
  AMREX_D_TERM(
    const Real fx = m_fac[0];, const Real fy = m_fac[1];
    , const Real fz = m_fac[2];);
  const Real detJ = AMREX_D_TERM(fx, *fy, *fz);

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    const auto& fac_ma = m_fac_fc[lev][d].arrays();
    const auto& det_ma = m_detJ_fc[lev][d].arrays();
    amrex::ParallelFor(
      m_fac_fc[lev][d], m_fac_fc[lev][d].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        auto f = fac_ma[box_no];
        auto dd = det_ma[box_no];
        AMREX_D_TERM(f(i, j, k, 0) = fx;, f(i, j, k, 1) = fy;
                     , f(i, j, k, 2) = fz;);
        dd(i, j, k) = detJ;
      });
  }
}
