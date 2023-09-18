#ifndef AMREX_MLABECCECLAP_1D_K_H
#define AMREX_MLABECCECLAP_1D_K_H

namespace amrex {

constexpr Real eps = 1.0e-6;

inline Real
getFaceStateUpwind(
  Real const& edgeVel, Real const& leftState, Real const& rightState)
{
  Real edstate = 0.0;
  if (edgeVel > eps) {
    edstate = leftState;
  }
  if (edgeVel < -eps) {
    edstate = rightState;
  }
  if (fabs(edgeVel) <= eps) {
    edstate = 0.5 * (leftState + rightState);
  }
  return edstate;
}

inline Real
getCellUpwindCoeff(Real const& leftedgeVel, Real const& rightedgeVel)
{

  Real coeff = 0.0;
  if (leftedgeVel * rightedgeVel < 0.0) { // Case where velocity change sign,
                                          // use centered + art. visc
    Real maxvel = std::max(fabs(leftedgeVel), fabs(rightedgeVel));
    coeff += 0.5 * (rightedgeVel - leftedgeVel); // centered part
    coeff += maxvel;                             // art. visc.
  } else { // same sign velocity, use upwind
    if (leftedgeVel < -eps)
      coeff -= leftedgeVel;
    else if (fabs(leftedgeVel) <= eps)
      coeff -= 0.5 * leftedgeVel;
    if (rightedgeVel > eps)
      coeff += rightedgeVel;
    else if (fabs(rightedgeVel) <= eps)
      coeff += 0.5 * rightedgeVel;
  }

  return coeff;
}

inline Real
getLeftCellUpwindCoeff(Real const& leftedgeVel, Real const& rightedgeVel)
{

  Real coeff = 0.0;
  if (leftedgeVel * rightedgeVel < 0.0) {
    Real maxvel = std::max(fabs(leftedgeVel), fabs(rightedgeVel));
    coeff -= 0.5 * leftedgeVel;
    coeff -= 0.5 * maxvel;
  } else {
    if (leftedgeVel > eps)
      coeff -= leftedgeVel;
    else if (fabs(leftedgeVel) < eps)
      coeff -= 0.5 * leftedgeVel;
  }
  return coeff;
}

inline Real
getRightCellUpwindCoeff(Real const& leftedgeVel, Real const& rightedgeVel)
{

  Real coeff = 0.0;
  if (leftedgeVel * rightedgeVel < 0.0) {
    Real maxvel = std::max(fabs(leftedgeVel), fabs(rightedgeVel));
    coeff += 0.5 * rightedgeVel;
    coeff += 0.5 * maxvel;
  } else {
    if (rightedgeVel < -eps)
      coeff += rightedgeVel;
    else if (fabs(rightedgeVel) < eps)
      coeff += 0.5 * rightedgeVel;
  }
  return coeff;
}

inline Real
getCellUpwindCoeff(
  Real const& leftedgeVel,
  Real const& rightedgeVel,
  Real const& relaxleft,
  Real const& relaxright)
{

  Real coeff = 0.0;
  if (
    leftedgeVel * relaxleft * rightedgeVel * relaxright <
    0.0) { // Case where velocity change sign, use centered + art. visc
    Real maxvel =
      std::max(fabs(leftedgeVel * relaxleft), fabs(rightedgeVel * relaxright));
    coeff += 0.5 * (rightedgeVel * relaxright -
                    leftedgeVel * relaxleft); // centered part
    coeff += maxvel;                          // art. visc.
  } else {                                    // same sign velocity, use upwind
    if (leftedgeVel < -eps)
      coeff -= leftedgeVel * relaxleft;
    else if (fabs(leftedgeVel) < eps)
      coeff -= 0.5 * leftedgeVel * relaxleft;
    if (rightedgeVel > eps)
      coeff += rightedgeVel * relaxright;
    else if (fabs(rightedgeVel) < eps)
      coeff += 0.5 * rightedgeVel * relaxright;
  }

  return coeff;
}

inline Real
getFaceUpwindFlux(
  Real const& edgeVel, Real const& leftState, Real const& rightState)
{

  Real Flux = 0.0;
  if (edgeVel > eps)
    Flux = edgeVel * leftState;
  else if (edgeVel < -eps)
    Flux = edgeVel * rightState;
  else
    Flux = edgeVel * 0.5 * (leftState + rightState);
  return Flux;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlabecceclap_adotx(
  Box const& box,
  Array4<Real> const& y,
  Array4<Real const> const& x,
  Array4<Real const> const& a,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
  Real alpha,
  Real beta,
  Real eta,
  int ncomp) noexcept
{
  const Real dhxsq = beta * dxinv[0] * dxinv[0];
  const Real dhx = eta * dxinv[0];

  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < ncomp; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      y(i, 0, 0, n) =
        alpha * a(i, 0, 0) * x(i, 0, 0, n) -
        dhxsq * (bX(i + 1, 0, 0) * (x(i + 1, 0, 0, n) - x(i, 0, 0, n)) -
                 bX(i, 0, 0) * (x(i, 0, 0, n) - x(i - 1, 0, 0, n))) -
        dhx * (cX(i + 1, 0, 0, n) *
                 getFaceStateUpwind(
                   cX(i + 1, 0, 0, n), x(i, 0, 0, n), x(i + 1, 0, 0, n)) -
               cX(i, 0, 0, n) *
                 getFaceStateUpwind(
                   cX(i, 0, 0, n), x(i - 1, 0, 0, n), x(i, 0, 0, n)));
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlabecceclap_adotx_os(
  Box const& box,
  Array4<Real> const& y,
  Array4<Real const> const& x,
  Array4<Real const> const& a,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  Array4<int const> const& osm,
  GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
  Real alpha,
  Real beta,
  Real eta,
  int ncomp) noexcept
{
  const Real dhxsq = beta * dxinv[0] * dxinv[0];
  const Real dhx = eta * dxinv[0];

  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < ncomp; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      if (osm(i, 0, 0) == 0) {
        y(i, 0, 0, n) = 0.0;
      } else {
        y(i, 0, 0, n) =
          alpha * a(i, 0, 0) * x(i, 0, 0, n) -
          dhxsq * (bX(i + 1, 0, 0) * (x(i + 1, 0, 0, n) - x(i, 0, 0, n)) -
                   bX(i, 0, 0) * (x(i, 0, 0, n) - x(i - 1, 0, 0, n))) -
          dhx * (cX(i + 1, 0, 0, n) *
                   getFaceStateUpwind(
                     cX(i + 1, 0, 0, n), x(i, 0, 0, n), x(i + 1, 0, 0, n)) -
                 cX(i, 0, 0, n) *
                   getFaceStateUpwind(
                     cX(i, 0, 0, n), x(i - 1, 0, 0, n), x(i, 0, 0, n)));
      }
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlabecceclap_normalize(
  Box const& box,
  Array4<Real> const& x,
  Array4<Real const> const& a,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  GpuArray<Real, AMREX_SPACEDIM> const& dxinv,
  Real alpha,
  Real beta,
  Real eta,
  int ncomp) noexcept
{
  const Real dhxsq = beta * dxinv[0] * dxinv[0];
  const Real dhx = eta * dxinv[0];

  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < ncomp; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      x(i, 0, 0, n) /=
        alpha * a(i, 0, 0) + dhxsq * (bX(i, 0, 0) + bX(i + 1, 0, 0)) -
        dhx * getCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n));
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlabecceclap_flux_x(
  Box const& box,
  Array4<Real> const& fx,
  Array4<Real const> const& sol,
  Array4<Real const> const& bx,
  Array4<Real const> const& cx,
  Real fac,
  Real facc,
  int ncomp) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < ncomp; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      fx(i, 0, 0, n) =
        -fac * bx(i, 0, 0) * (sol(i, 0, 0, n) - sol(i - 1, 0, 0, n)) -
        facc * getFaceUpwindFlux(
                 cx(i, 0, 0, n), sol(i - 1, 0, 0, n), sol(i, 0, 0, n));
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
mlabecceclap_flux_xface(
  Box const& box,
  Array4<Real> const& fx,
  Array4<Real const> const& sol,
  Array4<Real const> const& bx,
  Array4<Real const> const& cx,
  Real fac,
  Real facc,
  int xlen,
  int ncomp) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < ncomp; ++n) {
    int i = lo.x;
    fx(i, 0, 0, n) =
      -fac * bx(i, 0, 0) * (sol(i, 0, 0, n) - sol(i - 1, 0, 0, n)) -
      facc *
        getFaceUpwindFlux(cx(i, 0, 0, n), sol(i - 1, 0, 0, n), sol(i, 0, 0, n));
    i += xlen;
    fx(i, 0, 0, n) =
      -fac * bx(i, 0, 0) * (sol(i, 0, 0, n) - sol(i - 1, 0, 0, n)) -
      facc *
        getFaceUpwindFlux(cx(i, 0, 0, n), sol(i - 1, 0, 0, n), sol(i, 0, 0, n));
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
abeccec_gsrb(
  Box const& box,
  Array4<Real> const& phi,
  Array4<Real const> const& rhs,
  Real alpha,
  Real omega,
  Array4<Real const> const& a,
  Real dhxsq,
  Real dhx,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  Array4<int const> const& m0,
  Array4<int const> const& m1,
  Array4<Real const> const& f0,
  Array4<Real const> const& f1,
  Box const& vbox,
  int redblack,
  int nc) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
  const auto vlo = amrex::lbound(vbox);
  const auto vhi = amrex::ubound(vbox);

  for (int n = 0; n < nc; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      if ((i + redblack) % 2 == 0) {
        Real cf0 =
          (i == vlo.x and m0(vlo.x - 1, 0, 0) > 0) ? f0(vlo.x, 0, 0, n) : 0.0;
        Real cf1 =
          (i == vhi.x and m1(vhi.x + 1, 0, 0) > 0) ? f1(vhi.x, 0, 0, n) : 0.0;

        Real delta = dhxsq * (bX(i, 0, 0) * cf0 + bX(i + 1, 0, 0) * cf1);

        Real gamma =
          alpha * a(i, 0, 0) + dhxsq * (bX(i, 0, 0) + bX(i + 1, 0, 0)) -
          dhx * getCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n));

        Real rho =
          dhxsq * (bX(i, 0, 0) * phi(i - 1, 0, 0, n) +
                   bX(i + 1, 0, 0) * phi(i + 1, 0, 0, n)) +
          dhx * (getLeftCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n)) *
                   phi(i - 1, 0, 0, n) +
                 getRightCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n)) *
                   phi(i + 1, 0, 0, n));

        phi(i, 0, 0, n) = omega *
                            (rhs(i, 0, 0, n) + rho - phi(i, 0, 0, n) * delta) /
                            (gamma - delta) +
                          (1.0 - omega) * phi(i, 0, 0, n);
      }
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
abeccec_gsrb_os(
  Box const& box,
  Array4<Real> const& phi,
  Array4<Real const> const& rhs,
  Real alpha,
  Real omega,
  Array4<Real const> const& a,
  Real dhxsq,
  Real dhx,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  Array4<int const> const& m0,
  Array4<int const> const& m1,
  Array4<Real const> const& f0,
  Array4<Real const> const& f1,
  Array4<int const> const& osm,
  Box const& vbox,
  int redblack,
  int nc) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
  const auto vlo = amrex::lbound(vbox);
  const auto vhi = amrex::ubound(vbox);

  for (int n = 0; n < nc; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      if ((i + redblack) % 2 == 0) {
        if (osm(i, 0, 0) == 0) {
          phi(i, 0, 0) = 0.0;
        } else {
          Real cf0 =
            (i == vlo.x and m0(vlo.x - 1, 0, 0) > 0) ? f0(vlo.x, 0, 0, n) : 0.0;
          Real cf1 =
            (i == vhi.x and m1(vhi.x + 1, 0, 0) > 0) ? f1(vhi.x, 0, 0, n) : 0.0;

          Real delta = dhxsq * (bX(i, 0, 0) * cf0 + bX(i + 1, 0, 0) * cf1);

          Real gamma =
            alpha * a(i, 0, 0) + dhxsq * (bX(i, 0, 0) + bX(i + 1, 0, 0)) -
            dhx * getCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n));

          Real rho =
            dhxsq * (bX(i, 0, 0) * phi(i - 1, 0, 0, n) +
                     bX(i + 1, 0, 0) * phi(i + 1, 0, 0, n)) +
            dhx * (getLeftCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n)) *
                     phi(i - 1, 0, 0, n) +
                   getRightCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n)) *
                     phi(i + 1, 0, 0, n));

          phi(i, 0, 0, n) =
            omega * (rhs(i, 0, 0, n) + rho - phi(i, 0, 0, n) * delta) /
              (gamma - delta) +
            (1.0 - omega) * phi(i, 0, 0, n);
        }
      }
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
abeccec_gsrb_with_line_solve(
  Box const& box,
  Array4<Real> const& phi,
  Array4<Real const> const& rhs,
  Real alpha,
  Array4<Real const> const& a,
  Real dhx,
  Array4<Real const> const& bX,
  Array4<int const> const& m0,
  Array4<int const> const& m1,
  Array4<Real const> const& f0,
  Array4<Real const> const& f1,
  Box const& vbox,
  int redblack,
  int nc) noexcept
{
  amrex::Abort("abeccec_gsrb_with_line_solve not implemented in 1D");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE int
coarsen_overset_mask(
  Box const& bx,
  Array4<int> const& cmsk,
  Array4<int const> const& fmsk) noexcept
{
  int nerrors = 0;
  const auto lo = amrex::lbound(bx);
  const auto hi = amrex::ubound(bx);
  for (int i = lo.x; i <= hi.x; ++i) {
    int ii = 2 * i;
    cmsk(i, 0, 0) = fmsk(ii, 0, 0) + fmsk(ii + 1, 0, 0);
    if (cmsk(i, 0, 0) == 2) {
      cmsk(i, 0, 0) = 1;
    } else if (cmsk(i, 0, 0) != 0) {
      ++nerrors;
    }
  }
  return nerrors;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
overset_rescale_bcoef_x(
  Box const& box,
  Array4<Real> const& bX,
  Array4<int const> const& osm,
  int ncomp,
  Real osfac) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);
  for (int n = 0; n < ncomp; ++n) {
    for (int i = lo.x; i <= hi.x; ++i) {
      if ((osm(i - 1, 0, 0) + osm(i, 0, 0)) == 1) {
        bX(i, 0, 0, n) *= osfac;
      }
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
abeccec_chkdiag(
  Box const& box,
  Real alpha,
  Array4<Real const> const& a,
  Real dhxsq,
  Real dhx,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  int nc) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < nc; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      // Complete diag
      Real diag = alpha * a(i, 0, 0) -
                  dhxsq * (-bX(i, 0, 0, n) - bX(i + 1, 0, 0, n)) -
                  dhx * getCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n));

      // Offdiag
      Real offdiagdiff = -dhxsq * (bX(i, 0, 0, n) + bX(i + 1, 0, 0, n));

      Real offdiagadv =
        -dhx * (getLeftCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n)) +
                getRightCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n)));

      if (fabs(offdiagdiff + offdiagadv) > fabs(diag)) {
        // if ( i == 1 ) {
        amrex::Print() << " Abort LinOp not diagonal dominant at: " << i
                       << "\n";
        // amrex::Print() << j << " Val diag: " << diag << ", offdiagdiff/adv: "
        // << offdiagdiff << " " << offdiagadv << "\n"; amrex::Print() << j << "
        // Val diag: " << diag << ", left/right: " << left << "  " << right <<
        // "\n"; amrex::Print() << j << " " << fabs(diag) << " " <<
        // fabs(offdiagdiff+offdiagadv) << "\n"; amrex::Abort("MLABecCec is not
        // diagonal dominant !!!");
      }
    }
  }
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE void
abeccec_getdiag(
  Box const& box,
  Real alpha,
  Array4<Real const> const& a,
  Real dhxsq,
  Real dhx,
  Array4<Real const> const& bX,
  Array4<Real const> const& cX,
  int nc,
  Array4<Real> const& diag) noexcept
{
  const auto lo = amrex::lbound(box);
  const auto hi = amrex::ubound(box);

  for (int n = 0; n < nc; ++n) {
    AMREX_PRAGMA_SIMD
    for (int i = lo.x; i <= hi.x; ++i) {
      diag(i, 0, 0) =
        alpha * a(i, 0, 0) - dhxsq * (-bX(i, 0, 0, n) - bX(i + 1, 0, 0, n)) -
        dhx * getCellUpwindCoeff(cX(i, 0, 0, n), cX(i + 1, 0, 0, n));
    }
  }
}
} // namespace amrex
#endif
