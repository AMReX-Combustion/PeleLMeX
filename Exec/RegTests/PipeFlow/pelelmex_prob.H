#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMFData.H>
#include <PelePhysics.H>

struct MyProbParm : public ProbParmDefault
{
  amrex::Real T_mean = 298.0;
  amrex::Real meanFlowMag = 5.0;
  int meanFlowDir = 1;
  int flowType = 1;
  // Magnitude of the IC's transverse sin*sin perturbations.  The
  // historical default (5.0) deposits energy at wavenumbers up to
  // k=11 in x and k=4-5 in y/z, which is poorly resolved at modest
  // N (e.g. N=32) and pollutes self-convergence comparisons.  Set to
  // a smaller value (e.g. 0.0 for a clean uniform flow, or 0.5 for a
  // mild perturbation) when running a convergence study.
  amrex::Real perturbMag = 5.0;
  // Mesh-mapping scaling factors (read from ConstantMap.scaling_factor).
  // When mesh mapping is active the AMReX grid is Xi-space and physical
  // coordinates are  x_phys = prob_lo + (i+0.5) * dx * fac.  Defaults
  // of 1.0 reproduce the standard (unmapped) IC.
  amrex::Real fac_x = 1.0;
  amrex::Real fac_y = 1.0;
  amrex::Real fac_z = 1.0;
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
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

    // Convert AMReX indices to *physical* coordinates using the mesh-
    // mapping factors.  For the default fac = 1 this is identical to
    // the standard AMReX-coordinate IC.
    AMREX_D_TERM(
      const amrex::Real x =
        prob_lo[0] + (i + 0.5) * dx[0] * prob_parm.fac_x;
      , const amrex::Real y =
        prob_lo[1] + (j + 0.5) * dx[1] * prob_parm.fac_y;
      , const amrex::Real z =
        prob_lo[2] + (k + 0.5) * dx[2] * prob_parm.fac_z;);

    AMREX_D_TERM(
      const amrex::Real Lx = (prob_hi[0] - prob_lo[0]) * prob_parm.fac_x;
      , const amrex::Real Ly = (prob_hi[1] - prob_lo[1]) * prob_parm.fac_y;
      , const amrex::Real Lz = (prob_hi[2] - prob_lo[2]) * prob_parm.fac_z);

    auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    if (prob_parm.flowType == 1) {
      switch (prob_parm.meanFlowDir) {
      case 1:
        AMREX_D_TERM(
          state(i, j, k, VELX) = prob_parm.meanFlowMag;
          , state(i, j, k, VELY) = prob_parm.perturbMag *
                                   std::sin(2 * Pi * 4 * y / Ly) *
                                   std::sin(2 * Pi * 11 * x / Lx);
          , state(i, j, k, VELZ) = prob_parm.perturbMag *
                                   std::sin(2 * Pi * 5 * z / Lz) *
                                   std::sin(2 * Pi * 4 * x / Lx));
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
    } else if (prob_parm.flowType == 2) {
      // Turbulent-channel-like analytical profile (no EB needed):
      //   parabolic-ish mean profile in y across a channel of half-width Ly/2
      //   with superimposed divergence-free structured perturbations.
      const amrex::Real dpdx = 709.79; // Pa/m
      const amrex::Real c_f = 0.0062418;
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> u;
      u[0] = sqrt(dpdx * Ly / (0.4688 * c_f)) *
             (1.0 - pow(fabs(2.0 * y / Ly), 10.0));

      AMREX_D_TERM(
        state(i, j, k, VELX) =
          u[0] *
          (1.0 +
           std::sin(2 * Pi * 5 * y / Ly) * std::sin(2 * Pi * 14 * z / Lz) +
           std::sin(2 * Pi * 7 * y / Ly) *
             std::sin(2 * Pi * 6 * (z - 0.00541) / Lz));
        , state(i, j, k, VELY) =
            0.2 * u[0] *
            (std::sin(2 * Pi * 5 * y / Ly) *
               std::sin(2 * Pi * 7 * (x - 0.0013) / Lx) +
             std::sin(2 * Pi * 6 * y / Ly) * std::sin(2 * Pi * 13 * x / Lx));
        , state(i, j, k, VELZ) =
            0.2 * u[0] *
            (std::sin(2 * Pi * 11 * z / Lz) * std::sin(2 * Pi * 4 * x / Lx) +
             std::sin(2 * Pi * 3 * z / Lz) *
               std::sin(2 * Pi * 9 * (x - 0.00779) / Lx)));
    }

    if (!(is_incompressible == 0)) {
      return;
    }

    state(i, j, k, TEMP) = prob_parm.T_mean;

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

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int /*idir*/,
    const int /*sgn*/,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    auto eos = pele::physics::PhysicsType::eos();

    if (prob_parm.flowType == 1) {
      switch (prob_parm.meanFlowDir) {
      case 1:
        AMREX_D_TERM(
          s_ext[VELX] = prob_parm.meanFlowMag;, s_ext[VELY] = 0.0;
          , s_ext[VELZ] = 0.0);
        break;
      case 2:
        AMREX_D_TERM(
          s_ext[VELY] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
          , s_ext[VELZ] = 0.0);
        break;
      case 3:
        AMREX_D_TERM(
          s_ext[VELZ] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
          , s_ext[VELY] = 0.0);
        break;
      default:
        break;
      }

      massfrac[O2_ID] = 0.333;
      massfrac[N2_ID] = 0.667;

      s_ext[TEMP] = prob_parm.T_mean;

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] =
        RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
      }
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
