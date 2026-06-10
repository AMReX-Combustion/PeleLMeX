#include <PeleLMeX_TanhStretchMap.H>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>

using amrex::Real;

namespace {

// Below this magnitude the analytic tanh formulas are replaced by their
// identity limit (fac = 1, x_phys = xi, disp = 0).
constexpr Real BETA_EPS = Real(1.0e-8);

// fac(eta) along one axis with stretching parameter beta.  eta in
// [0, 1] is the normalized position; cells cluster at both ends when
// beta > 0 and the formula returns 1 when beta < BETA_EPS (identity).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
tanh_fac(Real eta, Real beta) noexcept
{
  if (std::abs(beta) < BETA_EPS) {
    return Real(1.0);
  }
  const Real arg = beta * (Real(2.0) * eta - Real(1.0));
  const Real sech = Real(1.0) / std::cosh(arg);
  return beta * sech * sech / std::tanh(beta);
}

// Offset from the low end of the axis, in physical-space units.  When
// scaled by 1/L it is the normalized position of x_phys; the actual
// x_phys is prob_lo + tanh_offset(eta, beta) * L.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
tanh_offset_norm(Real eta, Real beta) noexcept
{
  if (std::abs(beta) < BETA_EPS) {
    return eta;
  }
  return Real(0.5) *
         (Real(1.0) +
          std::tanh(beta * (Real(2.0) * eta - Real(1.0))) / std::tanh(beta));
}

} // namespace

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
  AMREX_ASSERT(lev >= 0 && lev < num_levels());

  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, Lxi{}, dxi{}, beta{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    Lxi[d] = geom.ProbHi(d) - geom.ProbLo(d);
    dxi[d] = geom.CellSize(d);
    beta[d] = m_beta[d];
  }

  // Evaluate fac on axis d at integer index n with offset off
  // (0 = lower node, 0.5 = cell center, 1 = upper node).
  auto fac_at = [=] AMREX_GPU_HOST_DEVICE(int d, int n, Real off) noexcept {
    const Real xi = plo[d] + (static_cast<Real>(n) + off) * dxi[d];
    const Real eta = (xi - plo[d]) / Lxi[d];
    return tanh_fac(eta, beta[d]);
  };

  // ---- Cell-centered: all axes at CC (offset 0.5) ----
  {
    const auto& fac_ma = m_fac_cc[lev].arrays();
    const auto& det_ma = m_detJ_cc[lev].arrays();
    amrex::ParallelFor(
      m_fac_cc[lev], m_fac_cc[lev].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        auto f = fac_ma[box_no];
        AMREX_D_TERM(
          const Real fx = fac_at(0, i, Real(0.5));
          , const Real fy = fac_at(1, j, Real(0.5));
          , const Real fz = fac_at(2, k, Real(0.5)););
        AMREX_D_TERM(f(i, j, k, 0) = fx;, f(i, j, k, 1) = fy;
                     , f(i, j, k, 2) = fz;);
        det_ma[box_no](i, j, k) = AMREX_D_TERM(fx, *fy, *fz);
      });
  }

  // ---- Nodal: all axes at ND (offset 0.0) ----
  {
    const auto& fac_ma = m_fac_nd[lev].arrays();
    const auto& det_ma = m_detJ_nd[lev].arrays();
    amrex::ParallelFor(
      m_fac_nd[lev], m_fac_nd[lev].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        auto f = fac_ma[box_no];
        AMREX_D_TERM(
          const Real fx = fac_at(0, i, Real(0.0));
          , const Real fy = fac_at(1, j, Real(0.0));
          , const Real fz = fac_at(2, k, Real(0.0)););
        AMREX_D_TERM(f(i, j, k, 0) = fx;, f(i, j, k, 1) = fy;
                     , f(i, j, k, 2) = fz;);
        det_ma[box_no](i, j, k) = AMREX_D_TERM(fx, *fy, *fz);
      });
  }

  // ---- Face-centered, per face direction d_face.  On axis d_face the
  // ---- position is at a node (offset 0); on the other axes it is cc
  // ---- (offset 0.5). ----
  for (int dface = 0; dface < AMREX_SPACEDIM; ++dface) {
    const int dface_cap = dface;
    const auto& fac_ma = m_fac_fc[lev][dface].arrays();
    const auto& det_ma = m_detJ_fc[lev][dface].arrays();
    amrex::ParallelFor(
      m_fac_fc[lev][dface], m_fac_fc[lev][dface].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        auto f = fac_ma[box_no];
        const Real off_x = (dface_cap == 0) ? Real(0.0) : Real(0.5);
#if (AMREX_SPACEDIM >= 2)
        const Real off_y = (dface_cap == 1) ? Real(0.0) : Real(0.5);
#endif
#if (AMREX_SPACEDIM == 3)
        const Real off_z = (dface_cap == 2) ? Real(0.0) : Real(0.5);
#endif
        AMREX_D_TERM(
          const Real fx = fac_at(0, i, off_x);
          , const Real fy = fac_at(1, j, off_y);
          , const Real fz = fac_at(2, k, off_z););
        AMREX_D_TERM(f(i, j, k, 0) = fx;, f(i, j, k, 1) = fy;
                     , f(i, j, k, 2) = fz;);
        det_ma[box_no](i, j, k) = AMREX_D_TERM(fx, *fy, *fz);
      });
  }

  amrex::Gpu::streamSynchronize();
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
  AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
  AMREX_ASSERT(disp_nd.nComp() >= AMREX_SPACEDIM);
  AMREX_ASSERT_WITH_MESSAGE(
    disp_nd.boxArray().ixType().nodeCentered(),
    "TanhStretchMap::fill_nodal_displacement: output MultiFab must be nodal");

  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, Lxi{}, dxi{}, beta{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    Lxi[d] = geom.ProbHi(d) - geom.ProbLo(d);
    dxi[d] = geom.CellSize(d);
    beta[d] = m_beta[d];
  }

  // disp on axis d at node index n: x_phys(node) - xi(node).
  auto disp_at = [=] AMREX_GPU_HOST_DEVICE(int d, int n) noexcept {
    const Real xi = plo[d] + static_cast<Real>(n) * dxi[d];
    const Real eta = (xi - plo[d]) / Lxi[d];
    const Real x_phys = plo[d] + Lxi[d] * tanh_offset_norm(eta, beta[d]);
    return x_phys - xi;
  };

  const auto& disp_ma = disp_nd.arrays();
  amrex::ParallelFor(
    disp_nd, amrex::IntVect(0),
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      auto d = disp_ma[box_no];
      AMREX_D_TERM(d(i, j, k, 0) = disp_at(0, i);
                   , d(i, j, k, 1) = disp_at(1, j);
                   , d(i, j, k, 2) = disp_at(2, k););
    });
  amrex::Gpu::streamSynchronize();
}
