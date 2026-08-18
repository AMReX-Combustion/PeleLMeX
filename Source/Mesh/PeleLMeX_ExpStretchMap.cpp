#include <PeleLMeX_ExpStretchMap.H>

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_Gpu.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

#include <string>

using amrex::Real;

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
  // L_xi == L_phys (length-preserving): the AMReX domain is the physical
  // domain and the mapping just redistributes cells within it.  The
  // non-stretched axes carry fac = 1.
  const int idir = m_dir;
  const int wall_end = m_wall_end;
  const Real beta = m_beta;
  const Real xi_lo = geom.ProbLo(idir);
  const Real L_xi = geom.ProbHi(idir) - xi_lo;

  fill_metrics(lev, geom, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
    return (d != idir)
             ? Real(1.0)
             : MeshMapEvaluator::exp_fac((xi - xi_lo) / L_xi, beta, wall_end);
  });
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
  const int idir = m_dir;
  const int wall_end = m_wall_end;
  const Real beta = m_beta;
  const Real xi_lo = geom.ProbLo(idir);
  const Real L_xi = geom.ProbHi(idir) - xi_lo;

  fill_nodal_disp_from(
    lev, geom, disp_nd, [=] AMREX_GPU_HOST_DEVICE(int d, Real xi) noexcept {
      return (d != idir)
               ? xi
               : xi_lo + L_xi * MeshMapEvaluator::exp_offset_norm(
                                  (xi - xi_lo) / L_xi, beta, wall_end);
    });
}
