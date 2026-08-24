#include <PeleLMeX_InteriorStretchMap.H>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <string>

using amrex::Real;

namespace {

// Below this magnitude a beta is snapped to exactly zero (uniform spacing
// on that side); with both sides zero the map is the identity.  The
// analytic formulas themselves stay well defined as beta -> 0 --- see
// MeshMapEvaluator::interior_sinhc --- so this is for reporting
// consistency with the other maps rather than numerical necessity.
constexpr Real BETA_EPS = Real(1.0e-8);

// Keep the split strictly inside (0,1) so the per-side normalisations
// s_c and 1-s_c never divide by zero.
constexpr Real SPLIT_EPS = Real(1.0e-12);

// Worst neighbouring-cell size ratio on a side carrying n_side cells,
// from d(ln fac)/ds = 2 beta tanh(beta eta) / W_side evaluated at the
// wall.  Asymptotic in n_side; accurate to <1% for n_side >~ 8.
Real
growth_estimate(Real beta, Real n_side) noexcept
{
  if (beta < BETA_EPS || n_side <= Real(0.0)) {
    return Real(1.0);
  }
  return std::exp(Real(2.0) * beta * std::tanh(beta) / n_side);
}

} // namespace

InteriorStretchMap::InteriorStretchMap()
{
  // ---- Domain, needed to turn a physical clustering point into a
  // ---- fraction of the axis.  geometry.prob_lo / prob_hi are required
  // ---- inputs for any AMReX run, and are identical on every level, so
  // ---- the mapping is fixed once here rather than per create_map() call.
  {
    amrex::ParmParse ppg("geometry");
    amrex::Vector<Real> plo(AMREX_SPACEDIM, 0.0);
    amrex::Vector<Real> phi(AMREX_SPACEDIM, 1.0);
    ppg.getarr("prob_lo", plo, 0, AMREX_SPACEDIM);
    ppg.getarr("prob_hi", phi, 0, AMREX_SPACEDIM);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      m_prob_lo[d] = plo[d];
      m_Lxi[d] = phi[d] - plo[d];
      if (m_Lxi[d] <= Real(0.0)) {
        amrex::Abort(
          "InteriorStretchMap: geometry.prob_hi must exceed geometry.prob_lo "
          "on axis " +
          std::to_string(d) + ".");
      }
    }
  }

  amrex::ParmParse pp("InteriorStretchMap");

  // ---- Coarsening rates ----
  amrex::Vector<Real> blo(AMREX_SPACEDIM, 0.0);
  amrex::Vector<Real> bhi(AMREX_SPACEDIM, 0.0);
  pp.queryarr("beta_lo", blo, 0, AMREX_SPACEDIM);
  pp.queryarr("beta_hi", bhi, 0, AMREX_SPACEDIM);
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (blo[d] < Real(0.0) || bhi[d] < Real(0.0)) {
      amrex::Abort(
        "InteriorStretchMap: beta_lo/beta_hi must be non-negative (got " +
        std::to_string(blo[d]) + " / " + std::to_string(bhi[d]) + " on axis " +
        std::to_string(d) + ").");
    }
    m_beta_lo[d] = (blo[d] < BETA_EPS) ? Real(0.0) : blo[d];
    m_beta_hi[d] = (bhi[d] < BETA_EPS) ? Real(0.0) : bhi[d];
  }

  // ---- Clustering point: physical coordinate wins over fraction ----
  amrex::Vector<Real> cfrac(AMREX_SPACEDIM, 0.5);
  pp.queryarr("center_frac", cfrac, 0, AMREX_SPACEDIM);
  if (pp.contains("center")) {
    amrex::Vector<Real> cphys(AMREX_SPACEDIM, 0.0);
    pp.getarr("center", cphys, 0, AMREX_SPACEDIM);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      cfrac[d] = (cphys[d] - m_prob_lo[d]) / m_Lxi[d];
    }
  }

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    const bool stretched =
      (m_beta_lo[d] > Real(0.0)) || (m_beta_hi[d] > Real(0.0));
    if (stretched && (cfrac[d] <= Real(0.0) || cfrac[d] >= Real(1.0))) {
      amrex::Abort(
        "InteriorStretchMap: the clustering point must lie strictly inside "
        "the domain on axis " +
        std::to_string(d) + " (got a fraction of " + std::to_string(cfrac[d]) +
        ").  Use ExpStretchMap for clustering at a boundary.");
    }
    m_center_frac[d] =
      std::fmin(std::fmax(cfrac[d], SPLIT_EPS), Real(1.0) - SPLIT_EPS);

    // Closed-form split, then the normalisation that follows from it.
    m_s_c[d] = std::fmin(
      std::fmax(
        MeshMapEvaluator::interior_split(
          m_center_frac[d], m_beta_lo[d], m_beta_hi[d]),
        SPLIT_EPS),
      Real(1.0) - SPLIT_EPS);
    m_N[d] = MeshMapEvaluator::interior_N(m_beta_lo[d], m_beta_hi[d], m_s_c[d]);
  }

  // ---- Start-up report, plus the one failure mode worth warning about:
  // ---- a side starved of cells has brutal cell-to-cell growth.
  Real max_growth = Real(1.15);
  pp.query("max_growth", max_growth);

  amrex::Vector<int> ncell(AMREX_SPACEDIM, 0);
  {
    amrex::ParmParse ppa("amr");
    ppa.queryarr("n_cell", ncell, 0, AMREX_SPACEDIM);
  }

  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (m_beta_lo[d] == Real(0.0) && m_beta_hi[d] == Real(0.0)) {
      continue;
    }
    const Real n = static_cast<Real>(ncell[d]);
    const Real n_lo = m_s_c[d] * n;
    const Real n_hi = (Real(1.0) - m_s_c[d]) * n;

    amrex::Print() << "   InteriorStretchMap axis " << d << ": beta = ("
                   << m_beta_lo[d] << ", " << m_beta_hi[d]
                   << "), x_c = " << center(d) << ", Xi_c = " << centerXi(d)
                   << " (s_c = " << m_s_c[d] << ")\n";
    if (ncell[d] > 0) {
      amrex::Print() << "                          cells " << int(n_lo + 0.5)
                     << " below / " << static_cast<int>(std::lround(n_hi))
                     << " above x_c"
                     << ", dx_min = " << m_Lxi[d] / (n * m_N[d])
                     << ", dx_max/dx_min = "
                     << std::fmax(
                          std::cosh(m_beta_lo[d]) * std::cosh(m_beta_lo[d]),
                          std::cosh(m_beta_hi[d]) * std::cosh(m_beta_hi[d]))
                     << "\n";

      const Real g = std::fmax(
        growth_estimate(m_beta_lo[d], n_lo),
        growth_estimate(m_beta_hi[d], n_hi));
      const bool starved = (n_lo < Real(4.0) && m_beta_lo[d] > Real(0.0)) ||
                           (n_hi < Real(4.0) && m_beta_hi[d] > Real(0.0));
      if (g > max_growth || starved) {
        amrex::Print()
          << "   WARNING: InteriorStretchMap axis " << d
          << " has a starved side: estimated worst neighbouring-cell size "
          << "ratio is " << g << " (limit " << max_growth << ").\n"
          << "            Reduce the larger beta, move the clustering point "
          << "towards it, or raise amr.n_cell on this axis;\n"
          << "            a side needs roughly 2*beta*tanh(beta)/ln(g_max) "
          << "cells to stay under g_max.\n";
      }
    }
  }
}

void
InteriorStretchMap::create_map(int lev, const amrex::Geometry& geom)
{
  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, Lxi{}, blo{}, bhi{}, sc{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    Lxi[d] = geom.ProbHi(d) - geom.ProbLo(d);
    blo[d] = m_beta_lo[d];
    bhi[d] = m_beta_hi[d];
    sc[d] = m_s_c[d];
  }

  fill_metrics(lev, geom, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
    return MeshMapEvaluator::interior_fac(
      blo[d], bhi[d], sc[d], (xi - plo[d]) / Lxi[d]);
  });
}

MeshMapEvaluator
InteriorStretchMap::make_evaluator() const
{
  MeshMapEvaluator e;
  e.m_kind = MeshMapEvaluator::Kind::InteriorStretch;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    e.m_p[d] = m_beta_lo[d];
    e.m_p2[d] = m_beta_hi[d];
    e.m_p3[d] = m_s_c[d];
    e.m_q[d] = -1;
  }
  return e;
}

void
InteriorStretchMap::fill_nodal_displacement(
  int lev, const amrex::Geometry& geom, amrex::MultiFab& disp_nd) const
{
  amrex::GpuArray<Real, AMREX_SPACEDIM> plo{}, Lxi{}, blo{}, bhi{}, sc{};
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    plo[d] = geom.ProbLo(d);
    Lxi[d] = geom.ProbHi(d) - geom.ProbLo(d);
    blo[d] = m_beta_lo[d];
    bhi[d] = m_beta_hi[d];
    sc[d] = m_s_c[d];
  }

  fill_nodal_disp_from(
    lev, geom, disp_nd, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
      return plo[d] + Lxi[d] * MeshMapEvaluator::interior_offset_norm(
                                 blo[d], bhi[d], sc[d], (xi - plo[d]) / Lxi[d]);
    });
}
