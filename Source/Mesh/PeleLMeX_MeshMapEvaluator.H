#ifndef PELELMEX_MESHMAP_EVALUATOR_H
#define PELELMEX_MESHMAP_EVALUATOR_H

#include <AMReX.H>
#include <AMReX_BLassert.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_SPACE.H>

#include <cmath>
#include <limits>

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
    TanhStretch = 3,
    InteriorStretch = 4
  };

  // Below this magnitude all stretching parameters revert to the
  // identity mapping (keeps the formulas numerically well-defined as
  // beta -> 0).
  static constexpr amrex::Real BETA_EPS = amrex::Real(1.0e-8);

  /// Mapping kind currently active.
  Kind m_kind{Kind::Identity};

  /// Per-axis primary parameter.  Meaning depends on \c m_kind:
  ///   Identity        : unused
  ///   Constant        : per-axis stretch factor
  ///   ExpStretch      : beta on the stretched axis, 0 elsewhere
  ///   TanhStretch     : per-axis beta (length-preserving)
  ///   InteriorStretch : per-axis beta_lo (coarsening rate towards prob_lo)
  amrex::Real m_p[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};

  /// Per-axis secondary Real parameter.  Meaning depends on \c m_kind:
  ///   InteriorStretch : per-axis beta_hi (coarsening rate towards prob_hi)
  ///   others          : unused (0)
  amrex::Real m_p2[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};

  /// Per-axis tertiary Real parameter.  Meaning depends on \c m_kind:
  ///   InteriorStretch : per-axis Xi-space split s_c in (0,1), i.e. the
  ///                     fraction of the cells placed below the clustering
  ///                     point.  Everything else about the map follows from
  ///                     (beta_lo, beta_hi, s_c) in closed form.
  ///   others          : unused (0)
  amrex::Real m_p3[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};

  /// Per-axis discrete parameter.  Meaning depends on \c m_kind:
  ///   ExpStretch  : 0 = wall at lo end, 1 = wall at hi end,
  ///                 -1 if this axis is not the stretched one
  ///   others      : unused (-1)
  int m_q[AMREX_SPACEDIM] = {AMREX_D_DECL(-1, -1, -1)};

  // ======================================================================
  //  Mapping primitives --- the single source of truth
  // ======================================================================
  //
  // Every concrete map's create_map() and fill_nodal_displacement(), and
  // every case of x_phys_from_xi() below, are expressed in terms of these.
  // They live here, in one place, precisely so that the metric a map writes
  // into its MultiFabs cannot drift away from the coordinate this evaluator
  // reports for the same point: for each family, *_fac is by construction
  // the derivative of *_offset_norm.
  //
  // Convention throughout: `s` is the Xi-space fraction of the axis,
  //     s = (xi - prob_lo) / L,  s in [0,1],
  // and *_offset_norm returns the physical fraction (x_phys - prob_lo)/L,
  // so that x_phys = prob_lo + L * offset_norm(s, ...).  All are
  // device-callable and well behaved as beta -> 0, where they reduce to
  // fac = 1 and offset_norm(s) = s (the identity map).

  // ---- ExpStretchMap primitives (one-sided exponential) ----------------

  /**
   * \brief \f$\mathrm{fac} = \beta e^{\beta\eta}/(e^\beta - 1)\f$, where
   *        \f$\eta\f$ measures distance from the wall end.
   *
   * \p wall_end is 0 for a wall at \c prob_lo, 1 for a wall at \c prob_hi
   * (the mirrored case).
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  exp_fac(amrex::Real s, amrex::Real beta, int wall_end) noexcept
  {
    if (std::abs(beta) < BETA_EPS) {
      return amrex::Real(1.0);
    }
    const amrex::Real eta = (wall_end == 1) ? (amrex::Real(1.0) - s) : s;
    return beta * std::exp(beta * eta) / (std::exp(beta) - amrex::Real(1.0));
  }

  /**
   * \brief \f$(x_{\mathrm{phys}} - x_{\mathrm{lo}})/L\f$ for the
   *        one-sided exponential map; the integral of \c exp_fac.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  exp_offset_norm(amrex::Real s, amrex::Real beta, int wall_end) noexcept
  {
    if (std::abs(beta) < BETA_EPS) {
      return s;
    }
    const amrex::Real eta = (wall_end == 1) ? (amrex::Real(1.0) - s) : s;
    const amrex::Real off = (std::exp(beta * eta) - amrex::Real(1.0)) /
                            (std::exp(beta) - amrex::Real(1.0));
    return (wall_end == 1) ? (amrex::Real(1.0) - off) : off;
  }

  // ---- TanhStretchMap primitives (two-sided, clusters at both ends) ----

  /**
   * \brief \f$\mathrm{fac} = \beta\,\mathrm{sech}^2(\beta(2s-1))/
   *        \tanh\beta\f$.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  tanh_fac(amrex::Real s, amrex::Real beta) noexcept
  {
    if (std::abs(beta) < BETA_EPS) {
      return amrex::Real(1.0);
    }
    const amrex::Real sech =
      amrex::Real(1.0) /
      std::cosh(beta * (amrex::Real(2.0) * s - amrex::Real(1.0)));
    return beta * sech * sech / std::tanh(beta);
  }

  /**
   * \brief \f$(x_{\mathrm{phys}} - x_{\mathrm{lo}})/L\f$ for the two-sided
   *        \c tanh map; the integral of \c tanh_fac.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  tanh_offset_norm(amrex::Real s, amrex::Real beta) noexcept
  {
    if (std::abs(beta) < BETA_EPS) {
      return s;
    }
    return amrex::Real(0.5) *
           (amrex::Real(1.0) +
            std::tanh(beta * (amrex::Real(2.0) * s - amrex::Real(1.0))) /
              std::tanh(beta));
  }

  // ---- InteriorStretchMap primitives -----------------------------------

  /** \brief \f$\sinh(z)/z\f$, accurate (and finite) as \f$z\to 0\f$. */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  interior_sinhc(amrex::Real z) noexcept
  {
    const amrex::Real z2 = z * z;
    if (z2 < amrex::Real(1.0e-8)) {
      return amrex::Real(1.0) +
             z2 * (amrex::Real(1.0) / amrex::Real(6.0) +
                   z2 * (amrex::Real(1.0) / amrex::Real(120.0)));
    }
    return std::sinh(z) / z;
  }

  /**
   * \brief \f$\mathcal{R}(\beta,\eta)=\int_0^\eta\cosh^2(\beta t)\,dt
   *        = \frac{\eta}{2} + \frac{\sinh(2\beta\eta)}{4\beta}\f$.
   *
   * \f$\mathcal{R}(\beta,\eta)\to\eta\f$ as \f$\beta\to 0\f$.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  interior_R(amrex::Real beta, amrex::Real eta) noexcept
  {
    return amrex::Real(0.5) * eta *
           (amrex::Real(1.0) + interior_sinhc(amrex::Real(2.0) * beta * eta));
  }

  /** \brief \f$\mathcal{G}(\beta)=\mathcal{R}(\beta,1)\f$. */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  interior_G(amrex::Real beta) noexcept
  {
    return interior_R(beta, amrex::Real(1.0));
  }

  /**
   * \brief Global normalisation \f$\mathcal{N}\f$ that makes the mapping
   *        preserve the total axis extent.  \c fac_min = 1/N.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  interior_N(amrex::Real beta_lo, amrex::Real beta_hi, amrex::Real s_c) noexcept
  {
    return s_c * interior_G(beta_lo) +
           (amrex::Real(1.0) - s_c) * interior_G(beta_hi);
  }

  /**
   * \brief Closed-form Xi-space split that puts the clustering point at the
   *        physical fraction \p p of the axis:
   * \f[
   *   s_c = \frac{p\,\mathcal{G}(\beta_{\mathrm{hi}})}
   *              {(1-p)\,\mathcal{G}(\beta_{\mathrm{lo}})
   *               + p\,\mathcal{G}(\beta_{\mathrm{hi}})}.
   * \f]
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real interior_split(
    amrex::Real p, amrex::Real beta_lo, amrex::Real beta_hi) noexcept
  {
    const amrex::Real G_lo = interior_G(beta_lo);
    const amrex::Real G_hi = interior_G(beta_hi);
    return p * G_hi / ((amrex::Real(1.0) - p) * G_lo + p * G_hi);
  }

  /**
   * \brief Stretch factor \f$dx_{\mathrm{phys}}/d\xi
   *        = \cosh^2(\beta_{\mathrm{side}}\eta)/\mathcal{N}\f$ at the
   *        Xi-space fraction \p s of the axis.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real interior_fac(
    amrex::Real beta_lo,
    amrex::Real beta_hi,
    amrex::Real s_c,
    amrex::Real s) noexcept
  {
    const amrex::Real a = (s <= s_c)
                            ? beta_lo * (s_c - s) / s_c
                            : beta_hi * (s - s_c) / (amrex::Real(1.0) - s_c);
    const amrex::Real c = std::cosh(a);
    return c * c / interior_N(beta_lo, beta_hi, s_c);
  }

  /**
   * \brief Physical-space fraction of the axis reached at the Xi-space
   *        fraction \p s, i.e. \f$(x_{\mathrm{phys}}-x_{\mathrm{lo}})/L\f$.
   *
   * Integrating \c interior_fac from 0 to \p s.  Returns \p s exactly when
   * both betas vanish.
   */
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real
  interior_offset_norm(
    amrex::Real beta_lo,
    amrex::Real beta_hi,
    amrex::Real s_c,
    amrex::Real s) noexcept
  {
    const amrex::Real N = interior_N(beta_lo, beta_hi, s_c);
    const amrex::Real xc_n = s_c * interior_G(beta_lo) / N;
    if (s <= s_c) {
      return xc_n - (s_c / N) * interior_R(beta_lo, (s_c - s) / s_c);
    }
    return xc_n + ((amrex::Real(1.0) - s_c) / N) *
                    interior_R(beta_hi, (s - s_c) / (amrex::Real(1.0) - s_c));
  }

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
      return plo + L * exp_offset_norm((xi - plo) / L, beta, m_q[idir]);
    }

    case Kind::TanhStretch: {
      const amrex::Real beta = m_p[idir];
      if (std::abs(beta) < BETA_EPS) {
        return xi;
      }
      const amrex::Real plo = gd.ProbLo()[idir];
      const amrex::Real L = gd.ProbHi()[idir] - plo;
      return plo + L * tanh_offset_norm((xi - plo) / L, beta);
    }

    case Kind::InteriorStretch: {
      const amrex::Real beta_lo = m_p[idir];
      const amrex::Real beta_hi = m_p2[idir];
      const amrex::Real s_c = m_p3[idir];
      if (std::abs(beta_lo) < BETA_EPS && std::abs(beta_hi) < BETA_EPS) {
        return xi;
      }
      const amrex::Real plo = gd.ProbLo()[idir];
      const amrex::Real L = gd.ProbHi()[idir] - plo;
      const amrex::Real s = (xi - plo) / L;
      return plo + L * interior_offset_norm(beta_lo, beta_hi, s_c, s);
    }
    }
    // Not reached: every Kind above returns.  The switch deliberately has
    // no default, so adding a Kind without handling it here is caught at
    // compile time by -Wswitch.  Should one slip through anyway, fail
    // loudly -- returning xi would silently degrade the run to an identity
    // map, which looks plausible in the output and is very hard to spot.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
      false, "MeshMapEvaluator::x_phys_from_xi: unhandled Kind -- a mesh map "
             "was added without a case in this switch");
    return std::numeric_limits<amrex::Real>::quiet_NaN();
  }
};

#endif // PELELMEX_MESHMAP_EVALUATOR_H
