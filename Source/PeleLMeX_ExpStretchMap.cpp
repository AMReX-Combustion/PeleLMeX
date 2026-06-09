#include <PeleLMeX_ExpStretchMap.H>

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

// Small-beta threshold below which the formulas are replaced by their
// identity limits (fac = 1, x_phys = xi, disp = 0).  This avoids the
// 0/0 at beta = 0 without chasing floating-point subtlety in user code.
constexpr Real BETA_EPS = Real(1.0e-8);

// fac(eta) = beta * exp(beta*eta) / (exp(beta) - 1), with eta in [0,1]
// measuring distance from the chosen wall end along the stretched axis.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
exp_fac(Real eta, Real beta) noexcept
{
  if (std::abs(beta) < BETA_EPS) {
    return Real(1.0);
  }
  const Real denom = std::exp(beta) - Real(1.0);
  return beta * std::exp(beta * eta) / denom;
}

// x_phys(eta) - xi_lo_wall, i.e. physical distance from the wall end
// along the stretched axis, where eta in [0,1] is the normalized
// distance from the wall.  Multiplied by L_xi gives the actual offset.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE Real
exp_offset(Real eta, Real beta) noexcept
{
  if (std::abs(beta) < BETA_EPS) {
    return eta;
  }
  const Real denom = std::exp(beta) - Real(1.0);
  return (std::exp(beta * eta) - Real(1.0)) / denom;
}

} // namespace

ExpStretchMap::ExpStretchMap()
{
  amrex::ParmParse pp("ExpStretchMap");
  pp.query("direction", m_dir);
  if (m_dir < 0 || m_dir >= AMREX_SPACEDIM) {
    amrex::Abort(
      "ExpStretchMap: 'direction' must be in [0, AMREX_SPACEDIM).  "
      "Set ExpStretchMap.direction to 0, 1, or 2.");
  }
  std::string wall{"lo"};
  pp.query("wall_end", wall);
  if (wall == "lo") {
    m_wall_end = 0;
  } else if (wall == "hi") {
    m_wall_end = 1;
  } else {
    amrex::Abort(
      "ExpStretchMap: 'wall_end' must be 'lo' or 'hi' (got '" + wall + "').");
  }
  pp.query("beta", m_beta);
  if (m_beta < Real(0.0)) {
    amrex::Abort(
      "ExpStretchMap: 'beta' must be non-negative (got " +
      std::to_string(m_beta) + ").");
  }
}

void
ExpStretchMap::create_map(int lev, const amrex::Geometry& geom)
{
  AMREX_ASSERT(lev >= 0 && lev < num_levels());

  // Geometry details for the stretched axis.  We use L_xi == L_phys
  // (length-preserving): the AMReX domain is the physical domain and
  // the mapping just redistributes cells within it.
  const int idir = m_dir;
  const int wall_end = m_wall_end;
  const Real beta = m_beta;
  const Real xi_lo = geom.ProbLo(idir);
  const Real xi_hi = geom.ProbHi(idir);
  const Real L_xi = xi_hi - xi_lo;
  const Real dxi = geom.CellSize(idir);

  // Helper: eta (normalized distance from wall) at an axis position xi.
  auto eta_of_xi = [=] AMREX_GPU_HOST_DEVICE(Real xi) noexcept -> Real {
    return (wall_end == 0) ? (xi - xi_lo) / L_xi : (xi_hi - xi) / L_xi;
  };

  // Helper: evaluate fac along the stretched axis at cell/node index n
  // with offset_in_cell in [0,1] (0 = node at low face of cell, 0.5 = cc,
  // 1 = node at high face).  For the non-stretched axes, fac = 1.
  auto fill_fac_array = [=] AMREX_GPU_HOST_DEVICE(
                          amrex::Array4<Real> const& fac, int i, int j, int k,
                          int stretched_index, Real offset_in_cell) noexcept {
    const Real xi_here =
      xi_lo + (static_cast<Real>(stretched_index) + offset_in_cell) * dxi;
    const Real eta = eta_of_xi(xi_here);
    const Real f = exp_fac(eta, beta);
    AMREX_D_TERM(fac(i, j, k, 0) = Real(1.0);, fac(i, j, k, 1) = Real(1.0);
                 , fac(i, j, k, 2) = Real(1.0););
    fac(i, j, k, idir) = f;
  };

  auto set_detJ = [=] AMREX_GPU_HOST_DEVICE(
                    amrex::Array4<Real> const& det, int i, int j, int k,
                    Real fac_stretched) noexcept {
    // detJ = product of fac_i; non-stretched directions contribute 1.
    det(i, j, k) = fac_stretched;
  };

  // ---- Cell-centered: index along stretched axis, offset 0.5 (cc). ----
  {
    const auto& fac_ma = m_fac_cc[lev].arrays();
    const auto& det_ma = m_detJ_cc[lev].arrays();
    amrex::ParallelFor(
      m_fac_cc[lev], m_fac_cc[lev].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        const int idx = (idir == 0) ? i : (idir == 1) ? j : k;
        fill_fac_array(fac_ma[box_no], i, j, k, idx, Real(0.5));
        const Real f = fac_ma[box_no](i, j, k, idir);
        set_detJ(det_ma[box_no], i, j, k, f);
      });
  }

  // ---- Nodal: index along stretched axis, offset 0.0 (node). ----
  {
    const auto& fac_ma = m_fac_nd[lev].arrays();
    const auto& det_ma = m_detJ_nd[lev].arrays();
    amrex::ParallelFor(
      m_fac_nd[lev], m_fac_nd[lev].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        const int idx = (idir == 0) ? i : (idir == 1) ? j : k;
        fill_fac_array(fac_ma[box_no], i, j, k, idx, Real(0.0));
        const Real f = fac_ma[box_no](i, j, k, idir);
        set_detJ(det_ma[box_no], i, j, k, f);
      });
  }

  // ---- Face-centered, per direction.  For face with normal = d,
  // the face position is at node (offset 0) along d and cc (offset 0.5)
  // along the other axes. ----
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    const auto& fac_ma = m_fac_fc[lev][d].arrays();
    const auto& det_ma = m_detJ_fc[lev][d].arrays();
    const int d_cap = d;
    amrex::ParallelFor(
      m_fac_fc[lev][d], m_fac_fc[lev][d].nGrowVect(),
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        // Offset along the stretched axis depends on whether the face
        // normal is the stretched axis (then offset = 0, at node) or
        // transverse (offset = 0.5, cc along the stretched axis).
        const int idx = (idir == 0) ? i : (idir == 1) ? j : k;
        const Real offset = (d_cap == idir) ? Real(0.0) : Real(0.5);
        fill_fac_array(fac_ma[box_no], i, j, k, idx, offset);
        const Real f = fac_ma[box_no](i, j, k, idir);
        set_detJ(det_ma[box_no], i, j, k, f);
      });
  }

  amrex::Gpu::streamSynchronize();
}

MeshMapEvaluator
ExpStretchMap::make_evaluator() const
{
  MeshMapEvaluator e;
  e.m_kind = MeshMapEvaluator::Kind::ExpStretch;
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    if (d == m_dir) {
      e.m_p[d] = m_beta;
      e.m_q[d] = m_wall_end; // 0 = lo, 1 = hi
    } else {
      e.m_p[d] = Real(0.0); // identity on the unstretched axes
      e.m_q[d] = -1;
    }
  }
  return e;
}

void
ExpStretchMap::fill_nodal_displacement(
  int lev, const amrex::Geometry& geom, amrex::MultiFab& disp_nd) const
{
  AMREX_ALWAYS_ASSERT(lev >= 0 && lev < num_levels());
  AMREX_ASSERT(disp_nd.nComp() >= AMREX_SPACEDIM);
  AMREX_ASSERT_WITH_MESSAGE(
    disp_nd.boxArray().ixType().nodeCentered(),
    "ExpStretchMap::fill_nodal_displacement: output MultiFab must be nodal");

  const int idir = m_dir;
  const int wall_end = m_wall_end;
  const Real beta = m_beta;
  const Real xi_lo = geom.ProbLo(idir);
  const Real xi_hi = geom.ProbHi(idir);
  const Real L_xi = xi_hi - xi_lo;
  const Real dxi = geom.CellSize(idir);

  const auto& disp_ma = disp_nd.arrays();
  amrex::ParallelFor(
    disp_nd, amrex::IntVect(0),
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      auto d = disp_ma[box_no];
      // Non-stretched directions: zero displacement.
      AMREX_D_TERM(d(i, j, k, 0) = Real(0.0);, d(i, j, k, 1) = Real(0.0);
                   , d(i, j, k, 2) = Real(0.0););

      // Stretched direction: evaluate the exact integral
      //   disp(xi) = x_phys(xi) - xi
      // at the node position.
      const int idx = (idir == 0) ? i : (idir == 1) ? j : k;
      const Real xi_here = xi_lo + static_cast<Real>(idx) * dxi;
      const Real eta =
        (wall_end == 0) ? (xi_here - xi_lo) / L_xi : (xi_hi - xi_here) / L_xi;
      // Offset (in physical space) from the wall along the stretched
      // axis: L_xi * (e^{beta*eta} - 1)/(e^beta - 1).  In the limit
      // beta -> 0 this reduces to L_xi * eta, i.e. x_phys == xi.
      const Real offset_from_wall = L_xi * exp_offset(eta, beta);
      const Real x_phys =
        (wall_end == 0) ? xi_lo + offset_from_wall : xi_hi - offset_from_wall;
      d(i, j, k, idir) = x_phys - xi_here;
    });
  amrex::Gpu::streamSynchronize();
}
