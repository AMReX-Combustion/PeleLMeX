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
ConstantMap::create_map(int lev, const amrex::Geometry& geom)
{
  amrex::GpuArray<Real, AMREX_SPACEDIM> fac{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    fac[d] = m_fac[d];
  }
  fill_metrics(
    lev, geom,
    [=] AMREX_GPU_HOST_DEVICE(int d, Real /*xi*/) noexcept { return fac[d]; });
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
ConstantMap::fill_nodal_displacement(
  int lev, const amrex::Geometry& geom, amrex::MultiFab& disp_nd) const
{
  // Linear map: x_phys = plo + fac * (xi - plo).  Equivalent to what the
  // former MeshMap base-class default computed, but stated explicitly
  // here rather than inherited by accident.
  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, fac{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    fac[d] = m_fac[d];
  }
  fill_nodal_disp_from(
    lev, geom, disp_nd, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
      return plo[d] + fac[d] * (xi - plo[d]);
    });
}
