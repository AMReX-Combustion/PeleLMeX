#ifndef PELELMEX_MESHMAP_EVALUATOR_H
#define PELELMEX_MESHMAP_EVALUATOR_H

#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#include <cmath>

/**
 * \brief Device-callable physical-coordinate evaluator for the mesh
 *        mapping.
 *
 * The mesh-mapping subsystem (see \c PeleLMeX_MeshMap.H) lays a
 * single-block, diagonal, axis-aligned mapping \f$X = \Phi(\xi)\f$ on
 * top of AMReX's uniform Xi-space grid.  Numerical operators consume
 * the metric MultiFabs (\c fac, \c detJ); user code in
 * \c initdata_mapped / \c bcnormal often wants the physical
 * coordinate itself.  This struct is the lightweight, GPU-clean
 * conduit between the two worlds.
 *
 * It is a trivially-copyable POD: an enum tag identifying which
 * concrete mapping is active plus a small fixed-size payload of
 * per-axis parameters.  The host-side factory
 * \c MeshMap::make_evaluator() returns a populated instance; the
 * driver caches one on \c PeleLM and passes it by value into device
 * lambdas.  No virtual dispatch, no MultiFab lookups per call.
 *
 * When no mesh mapping is active the evaluator defaults to
 * \c Kind::Identity and \c x_phys_cc / \c x_phys_face return the
 * unmapped \f$\xi\f$-space coordinate, i.e. the formula
 * \f$x = \text{prob\_lo} + (\text{idx}+0.5)\,\Delta x\f$.
 */
struct MeshMapEvaluator
{
  enum class Kind : int {
    Identity = 0,
    Constant = 1,
    ExpStretch = 2,
    TanhStretch = 3
  };

  // Below this magnitude all stretching parameters revert to the
  // identity mapping (keeps the formulas numerically well-defined as
  // beta -> 0).
  static constexpr amrex::Real BETA_EPS = amrex::Real(1.0e-8);

  /// Mapping kind currently active.
  Kind m_kind{Kind::Identity};

  /// Per-axis primary parameter.  Meaning depends on \c m_kind:
  ///   Identity    : unused
  ///   Constant    : per-axis stretch factor
  ///   ExpStretch  : beta on the stretched axis, 0 elsewhere
  ///   TanhStretch : per-axis beta (length-preserving)
  amrex::Real m_p[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};

  /// Per-axis discrete parameter.  Meaning depends on \c m_kind:
  ///   ExpStretch  : 0 = wall at lo end, 1 = wall at hi end,
  ///                 -1 if this axis is not the stretched one
  ///   others      : unused (-1)
  int m_q[AMREX_SPACEDIM] = {AMREX_D_DECL(-1, -1, -1)};

  /**
   * \brief Physical coordinate of the cell centre with integer index
   *        \p idx along axis \p idir.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  x_phys_cc(int idir, int idx, amrex::GeometryData const& gd) const noexcept
  {
    const amrex::Real xi =
      gd.ProbLo()[idir] +
      (static_cast<amrex::Real>(idx) + amrex::Real(0.5)) * gd.CellSize()[idir];
    return x_phys_from_xi(idir, xi, gd);
  }

  /**
   * \brief Physical size of domain along axis \p idir.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  domain_size_phys(int idir, amrex::GeometryData const& gd) const noexcept
  {
    return x_phys_from_xi(idir, gd.ProbHi()[idir], gd) -
           x_phys_from_xi(idir, gd.ProbLo()[idir], gd);
  }

  /**
   * \brief Physical coordinate of the face between cells \p idx-1 and
   *        \p idx along axis \p idir.
   *
   * \p idx is the nodal index (face position), so face 0 sits at
   * \c prob_lo and face N at \c prob_hi.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  x_phys_face(int idir, int idx, amrex::GeometryData const& gd) const noexcept
  {
    const amrex::Real xi =
      gd.ProbLo()[idir] + static_cast<amrex::Real>(idx) * gd.CellSize()[idir];
    return x_phys_from_xi(idir, xi, gd);
  }

  /**
   * \brief Lower-level helper: apply the mapping to an arbitrary
   *        \f$\xi\f$-space coordinate along axis \p idir.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real x_phys_from_xi(
    int idir, amrex::Real xi, amrex::GeometryData const& gd) const noexcept
  {
    switch (m_kind) {
    case Kind::Identity:
      return xi;

    case Kind::Constant: {
      // Linear: X = plo + fac * (xi - plo).
      const amrex::Real plo = gd.ProbLo()[idir];
      return plo + m_p[idir] * (xi - plo);
    }

    case Kind::ExpStretch: {
      const amrex::Real beta = m_p[idir];
      if (std::abs(beta) < BETA_EPS) {
        return xi;
      }
      const amrex::Real plo = gd.ProbLo()[idir];
      const amrex::Real L = gd.ProbHi()[idir] - plo;
      amrex::Real eta = (xi - plo) / L;
      // wall_end == 1 (wall at hi): mirror eta about 0.5 first,
      // apply the lo-wall formula, then mirror the result back.
      if (m_q[idir] == 1) {
        eta = amrex::Real(1.0) - eta;
      }
      const amrex::Real num = std::exp(beta * eta) - amrex::Real(1.0);
      const amrex::Real den = std::exp(beta) - amrex::Real(1.0);
      amrex::Real off = num / den;
      if (m_q[idir] == 1) {
        off = amrex::Real(1.0) - off;
      }
      return plo + L * off;
    }

    case Kind::TanhStretch: {
      const amrex::Real beta = m_p[idir];
      if (std::abs(beta) < BETA_EPS) {
        return xi;
      }
      const amrex::Real plo = gd.ProbLo()[idir];
      const amrex::Real L = gd.ProbHi()[idir] - plo;
      const amrex::Real eta = (xi - plo) / L;
      const amrex::Real off =
        amrex::Real(0.5) *
        (amrex::Real(1.0) +
         std::tanh(beta * (amrex::Real(2.0) * eta - amrex::Real(1.0))) /
           std::tanh(beta));
      return plo + L * off;
    }
    }
    // Unreachable; silences compiler warnings on some toolchains.
    return xi;
  }
};

#endif // PELELMEX_MESHMAP_EVALUATOR_H
