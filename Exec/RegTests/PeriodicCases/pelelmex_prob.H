#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PelePhysics.H>

struct MyProbParm : public ProbParmDefault
{
  int probType = 0;
  amrex::Real T_mean = 298.0;
  amrex::Real rvort = 0.07;
  amrex::Real xvort = 0.5;
  amrex::Real yvort = 0.5;
  amrex::Real forcevort = 6.0;
  int meanFlowDir = 1;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real rgauss = 0.1;
  amrex::Real xgauss = 0.5;
  amrex::Real ygauss = 0.5;
  amrex::Real ampgauss = 0.1;
  int gauss_type = 0;
  amrex::Real gaussTime = 0.001;
  amrex::Real gaussDiff = 0.000001;

  // ShearLayer: temporally-developing planar mixing layer.  Periodic in x,
  // bounded in y by slip walls placed far from the layer, so the only
  // feature in the domain is the interior shear layer at y = shear_y0 --
  // which is what InteriorStretchMap is for.
  amrex::Real shear_y0 = 0.0;    // physical y of the layer
  amrex::Real shear_dw = 0.004;  // vorticity thickness (du/dy|max = dU/dw)
  amrex::Real shear_dU = 1.0;    // velocity difference across the layer
  amrex::Real shear_Uc = 0.0;    // mean (convective) velocity
  amrex::Real shear_pert = 0.02; // max|v'| as a fraction of dU
  int shear_nwave = 1;           // KH wavelengths across the periodic x extent
  amrex::Real shear_ampY = 0.0;  // optional scalar step amplitude (low-Mach)
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  // Every problem type here uses a y coordinate, so this case has always
  // been 2D/3D only.  Say so, rather than reading prob_lo[1] out of bounds.
  static_assert(
    AMREX_SPACEDIM >= 2, "PeriodicCases requires AMREX_SPACEDIM >= 2");

  /**
   * \brief Set the state at (i,j,k) from the *physical* coordinates of the
   *        cell centre.
   *
   * Shared by initdata() and initdata_mapped() so the two entry points
   * cannot disagree: the only difference between them is how (x, y) and the
   * physical x-extent Lx are obtained.  On an unmapped run the mapped path
   * reduces to the uniform-grid formulas exactly (the evaluator falls back
   * to Kind::Identity), so this refactor is bit-neutral there.
   */
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void set_state(
    int i,
    int j,
    int k,
    int is_incompressible,
    amrex::Array4<amrex::Real> const& state,
    amrex::Real x,
    amrex::Real y,
    amrex::Real Lx,
    MyProbParm const& prob_parm)
  {
    amrex::ignore_unused(x, y, Lx);
    auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    // Species IDs
#ifdef O2_ID
    int sp1_ID = O2_ID;
#endif
#ifdef N2a_ID
    int sp1_ID = N2a_ID;
#endif
    int sp2_ID = N2_ID;

    massfrac[sp1_ID] = 0.233;
    massfrac[sp2_ID] = 0.767;

    if (prob_parm.probType == 0) { // CoVo
      const amrex::Real deltax = x - prob_parm.xvort;
      const amrex::Real deltay = y - prob_parm.yvort;
      const amrex::Real d_sq = deltax * deltax + deltay * deltay;
      const amrex::Real r_sq = prob_parm.rvort * prob_parm.rvort;
      const amrex::Real u_vort =
        -prob_parm.forcevort * deltay / r_sq * exp(-d_sq / r_sq / 2.);
      const amrex::Real v_vort =
        prob_parm.forcevort * deltax / r_sq * exp(-d_sq / r_sq / 2.);

      switch (prob_parm.meanFlowDir) {
      case 1:
        AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag + u_vort;
                     , state(i, j, k, VELY) = v_vort;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -1:
        AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag + u_vort;
                     , state(i, j, k, VELY) = v_vort;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case 2:
        AMREX_D_TERM(state(i, j, k, VELX) = u_vort;
                     , state(i, j, k, VELY) = prob_parm.meanFlowMag + v_vort;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -2:
        AMREX_D_TERM(state(i, j, k, VELX) = u_vort;
                     , state(i, j, k, VELY) = -prob_parm.meanFlowMag + v_vort;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case 3:
        AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag + u_vort;
                     , state(i, j, k, VELY) = prob_parm.meanFlowMag + v_vort;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -3:
        AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag + u_vort;
                     , state(i, j, k, VELY) = -prob_parm.meanFlowMag + v_vort;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      default:
        break;
      }

      if (is_incompressible != 0) {
        return;
      }

      state(i, j, k, TEMP) = prob_parm.T_mean;
    } else if (prob_parm.probType == 1) { // CoGau

      switch (prob_parm.meanFlowDir) {
      case 1:
        AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
        break;
      case -1:
        AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
        break;
      case 2:
        AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                     , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -2:
        AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                     , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case 3:
        AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -3:
        AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      default:
        break;
      }

      if (is_incompressible != 0) {
        return;
      }

      state(i, j, k, TEMP) = prob_parm.T_mean;
      const amrex::Real deltax = x - prob_parm.xgauss;
      const amrex::Real deltay = y - prob_parm.ygauss;
      const amrex::Real d_sq = deltax * deltax + deltay * deltay;
      const amrex::Real r_sq = prob_parm.rgauss * prob_parm.rgauss;
      if (prob_parm.gauss_type == 0) { // Spec
        massfrac[sp1_ID] +=
          massfrac[sp1_ID] * prob_parm.ampgauss * std::exp(-d_sq / r_sq);
        massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
      } else if (prob_parm.gauss_type == 1) { // Temp
        state(i, j, k, TEMP) =
          prob_parm.T_mean *
          (1.0 + prob_parm.ampgauss * std::exp(-d_sq / r_sq));
      } else {
        amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
      }
    } else if (prob_parm.probType == 3) { // CoTanH

      switch (prob_parm.meanFlowDir) {
      case 1:
        AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
        break;
      case -1:
        AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
        break;
      case 2:
        AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                     , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -2:
        AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                     , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case 3:
        AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      case -3:
        AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                     , state(i, j, k, VELZ) = 0.0);
        break;
      default:
        break;
      }

      if (is_incompressible != 0) {
        return;
      }

      state(i, j, k, TEMP) = prob_parm.T_mean;
      const amrex::Real deltax = x - prob_parm.xgauss;
      const amrex::Real deltay = y - prob_parm.ygauss;
      const amrex::Real rad = std::sqrt(deltax * deltax + deltay * deltay);
      const amrex::Real eta =
        0.5 *
        (1.0 -
         std::tanh(2.0 * (rad - prob_parm.rgauss) / (0.05 * prob_parm.rgauss)));
      if (prob_parm.gauss_type == 0) { // Spec
        massfrac[sp1_ID] += massfrac[sp1_ID] * prob_parm.ampgauss * eta;
        massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
      } else if (prob_parm.gauss_type == 1) { // Temp
        state(i, j, k, TEMP) =
          prob_parm.T_mean * (1.0 + prob_parm.ampgauss * eta);
      } else {
        amrex::Abort("Unknown tanh_type: should be either Temp or Spec");
      }
    } else if (prob_parm.probType == 2) { // DiffGau

      AMREX_D_TERM(state(i, j, k, VELX) = 0.0;, state(i, j, k, VELY) = 0.0;
                   , state(i, j, k, VELZ) = 0.0);

      if (is_incompressible != 0) {
        return;
      }

      state(i, j, k, TEMP) = prob_parm.T_mean;
      // Set up the Gaussian as the solution of the diffusion of a delta dirac
      // initial distribution after time gaussTime of a diffusion process with
      // constant diffusion coeff. gaussDiff. Use ampgauss to avoid
      // undershoot/overshoot of the 2-species mixture.
      const amrex::Real deltax = x - prob_parm.xgauss;
      amrex::Real denom_inv =
        prob_parm.ampgauss /
        std::sqrt(4 * Pi * prob_parm.gaussTime * prob_parm.gaussDiff);
      const amrex::Real d_sq = deltax * deltax;
      const amrex::Real r_sq = 4.0 * prob_parm.gaussTime * prob_parm.gaussDiff;
      if (prob_parm.gauss_type == 0) { // Spec
        massfrac[sp1_ID] = denom_inv * std::exp(-d_sq / r_sq);
        massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
      } else if (prob_parm.gauss_type == 1) { // Temp
        // TODO: find a better way, one with an analytical solution
        state(i, j, k, TEMP) =
          prob_parm.T_mean *
          (1.0 + prob_parm.ampgauss * std::exp(-d_sq / r_sq));
      } else {
        amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
      }
    } else if (prob_parm.probType == 4) { // ShearLayer

      // Base profile, with dw the vorticity thickness so that
      // du/dy|max = dU/dw at y = y0:
      //   u(y) = Uc + (dU/2) tanh(2 (y - y0) / dw)
      const amrex::Real dy = y - prob_parm.shear_y0;
      const amrex::Real th = std::tanh(2.0 * dy / prob_parm.shear_dw);
      amrex::Real u_shear = prob_parm.shear_Uc + 0.5 * prob_parm.shear_dU * th;
      amrex::Real v_shear = 0.0;

      // Kelvin-Helmholtz seed derived from a streamfunction, so the
      // perturbation is analytically divergence-free:
      //   psi' = A sin(kx x) exp(-(dy/dw)^2),  kx = 2 pi n / Lx
      //   u'   = -d psi'/dy,   v' = +d psi'/dx
      // A is scaled so that max|v'| = shear_pert * dU.
      if (prob_parm.shear_pert > 0.0 && prob_parm.shear_nwave > 0) {
        const amrex::Real kx =
          2.0 * Pi * static_cast<amrex::Real>(prob_parm.shear_nwave) / Lx;
        const amrex::Real amp = prob_parm.shear_pert * prob_parm.shear_dU / kx;
        const amrex::Real dw2 = prob_parm.shear_dw * prob_parm.shear_dw;
        const amrex::Real env = std::exp(-dy * dy / dw2);
        u_shear += amp * std::sin(kx * x) * (2.0 * dy / dw2) * env;
        v_shear += amp * kx * std::cos(kx * x) * env;
      }

      AMREX_D_TERM(state(i, j, k, VELX) = u_shear;
                   , state(i, j, k, VELY) = v_shear;
                   , state(i, j, k, VELZ) = 0.0);

      if (is_incompressible != 0) {
        return;
      }

      // Low-Mach: tag the two streams with a scalar step on the same profile
      // so the mixing is visible in a passive-ish scalar as well as in
      // vorticity.
      state(i, j, k, TEMP) = prob_parm.T_mean;
      if (prob_parm.shear_ampY != 0.0) {
        massfrac[sp1_ID] *= 1.0 + prob_parm.shear_ampY * 0.5 * (1.0 + th);
        massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
      }
    }

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
    state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
    }
  }

  /// Legacy entry point: uniform-grid cell-centre coordinates.
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata(
    int i,
    int j,
    int k,
    int is_incompressible,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();
    const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
    const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
    set_state(
      i, j, k, is_incompressible, state, x, y, prob_hi[0] - prob_lo[0],
      prob_parm);
  }

  /**
   * \brief Mapping-aware entry point (preferred by the driver when present).
   *
   * Evaluates the IC at the *physical* cell-centre position, which is what
   * makes the stretched-mesh runs meaningful: under InteriorStretchMap the
   * uniform-grid formula would place the shear layer at the wrong y.  Lx is
   * taken from the evaluator too, so the KH seed keeps an integer number of
   * wavelengths in the physical domain.
   */
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata_mapped(
    int i,
    int j,
    int k,
    int is_incompressible,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MeshMapEvaluator const& mmap,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real x = mmap.x_phys_cc(0, i, geomdata);
    const amrex::Real y = mmap.x_phys_cc(1, j, geomdata);
    set_state(
      i, j, k, is_incompressible, state, x, y,
      mmap.domain_size_phys(0, geomdata), prob_parm);
  }

  // Only used in non-periodic inflow cases
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    if (idir == (prob_parm.meanFlowDir - 1) && sgn == 1) {

      auto eos = pele::physics::PhysicsType::eos();
      // Species IDs
#ifdef O2_ID
      int sp1_ID = O2_ID;
#endif
#ifdef N2a_ID
      int sp1_ID = N2a_ID;
#endif
      int sp2_ID = N2_ID;

      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      massfrac[sp1_ID] = 0.233;
      massfrac[sp2_ID] = 0.767;

      s_ext[TEMP] = prob_parm.T_mean;
      amrex::Real P_cgs = prob_parm.P_mean * 10.0;
      // Density
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      // Enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(s_ext[TEMP], massfrac, h_cgs);
      s_ext[RHOH] = h_cgs * 1.0e-4 * s_ext[DENSITY];

      // Species mass
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
      }
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
