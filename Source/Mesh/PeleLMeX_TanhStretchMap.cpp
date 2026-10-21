#include <PeleLMeX_TanhStretchMap.H>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

using amrex::Real;

TanhStretchMap::TanhStretchMap()
{
  amrex::ParmParse pp("TanhStretchMap");
  amrex::Vector<amrex::Real> bvec(AMREX_SPACEDIM, 0.0);
  pp.queryarr("beta", bvec, 0, AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (bvec[d] < Real(0.0)) {
      amrex::Abort(
        "TanhStretchMap: beta values must be non-negative (got " +
        std::to_string(bvec[d]) + " on axis " + std::to_string(d) + ").");
    }
    m_beta[d] = bvec[d];
  }
}

void
TanhStretchMap::create_map(int lev, const amrex::Geometry& geom)
{
  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, Lxi{}, beta{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    Lxi[d] = geom.ProbHi(d) - geom.ProbLo(d);
    beta[d] = m_beta[d];
  }

  fill_metrics(lev, geom, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
    return MeshMapEvaluator::tanh_fac((xi - plo[d]) / Lxi[d], beta[d]);
  });
}

MeshMapEvaluator
TanhStretchMap::make_evaluator() const
{
  MeshMapEvaluator e;
  e.m_kind = MeshMapEvaluator::Kind::TanhStretch;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    e.m_p[d] = m_beta[d];
    e.m_q[d] = -1;
  }
  return e;
}

void
TanhStretchMap::fill_nodal_displacement(
  int lev, const amrex::Geometry& geom, amrex::MultiFab& disp_nd) const
{
  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, Lxi{}, beta{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    Lxi[d] = geom.ProbHi(d) - geom.ProbLo(d);
    beta[d] = m_beta[d];
  }

  fill_nodal_disp_from(
    lev, geom, disp_nd, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
      return plo[d] + Lxi[d] * MeshMapEvaluator::tanh_offset_norm(
                                 (xi - plo[d]) / Lxi[d], beta[d]);
    });
}
