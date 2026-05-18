#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMFData.H>
#include <PelePhysics.H>

struct MyProbParm : public ProbParmDefault
{
  amrex::Real T_mean = 298.0;
  amrex::Real meanFlowMag = 10.0;
  // Fractional amplitude of the initial-condition velocity perturbation
  // relative to meanFlowMag. 0.0 -> uniform IC. The perturbation is a
  // deterministic sum of sinusoids spanning a range of wavenumbers; the
  // backward-facing step's shear-layer instability scrambles it into
  // broadband turbulence within a flow-through or two.
  amrex::Real turb_seed_amplitude = 0.05;
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

    AMREX_D_TERM(
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
      , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

    AMREX_D_TERM(
      const amrex::Real Lx = prob_hi[0] - prob_lo[0];
      , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
      , const amrex::Real Lz = prob_hi[2] - prob_lo[2];);

    auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    // Only seed perturbation in the fluid region above/ahead of the step
    // (i.e., not inside the EB box); cells covered by EB get overwritten
    // later anyway, but keeping the seed clean is cheap and avoids weird
    // derived-plot artifacts near the cut cells.
    const bool in_step_block = (x < 0.01) && (y < 0.0);

    const amrex::Real U = prob_parm.meanFlowMag;
    const amrex::Real A = prob_parm.turb_seed_amplitude * U;

    // Multi-mode incompressible-looking velocity perturbation. Three
    // independent triads per component with incommensurate wavenumbers
    // and phase offsets so the pattern doesn't reduce to a simple
    // standing wave along any axis. Amplitude tapered to zero at the
    // wall-normal boundaries by (1 - (2y/Ly)^2).
    const amrex::Real ywall = 1.0 - std::pow(2.0 * y / Ly, 2.0);
    const amrex::Real wy = amrex::max<amrex::Real>(0.0, ywall);

    const amrex::Real pu =
      (std::sin(2.0 * Pi * 3.0 * x / Lx + 0.17) *
         std::cos(2.0 * Pi * 5.0 * y / Ly)
#if AMREX_SPACEDIM == 3
         * std::sin(2.0 * Pi * 4.0 * z / Lz)
#endif
       + std::sin(2.0 * Pi * 7.0 * x / Lx - 0.43) *
           std::cos(2.0 * Pi * 3.0 * y / Ly + 0.91)
#if AMREX_SPACEDIM == 3
           * std::sin(2.0 * Pi * 2.0 * z / Lz)
#endif
       + std::sin(2.0 * Pi * 11.0 * x / Lx + 1.23) *
           std::cos(2.0 * Pi * 2.0 * y / Ly)
#if AMREX_SPACEDIM == 3
           * std::sin(2.0 * Pi * 7.0 * z / Lz)
#endif
      );
    const amrex::Real pv =
      (std::cos(2.0 * Pi * 4.0 * x / Lx) *
         std::sin(2.0 * Pi * 4.0 * y / Ly + 0.31)
#if AMREX_SPACEDIM == 3
         * std::cos(2.0 * Pi * 3.0 * z / Lz + 0.77)
#endif
       + std::cos(2.0 * Pi * 6.0 * x / Lx + 0.55) *
           std::sin(2.0 * Pi * 2.0 * y / Ly)
#if AMREX_SPACEDIM == 3
           * std::cos(2.0 * Pi * 5.0 * z / Lz)
#endif
       + std::cos(2.0 * Pi * 9.0 * x / Lx) *
           std::sin(2.0 * Pi * 3.0 * y / Ly - 0.22)
#if AMREX_SPACEDIM == 3
           * std::cos(2.0 * Pi * 6.0 * z / Lz + 1.05)
#endif
      );
#if AMREX_SPACEDIM == 3
    const amrex::Real pw =
      (std::sin(2.0 * Pi * 5.0 * x / Lx + 0.65) *
         std::cos(2.0 * Pi * 3.0 * y / Ly) *
         std::sin(2.0 * Pi * 3.0 * z / Lz - 0.38) +
       std::sin(2.0 * Pi * 8.0 * x / Lx) *
         std::cos(2.0 * Pi * 4.0 * y / Ly + 0.47) *
         std::sin(2.0 * Pi * 4.0 * z / Lz) +
       std::sin(2.0 * Pi * 2.0 * x / Lx - 0.93) *
         std::cos(2.0 * Pi * 5.0 * y / Ly) *
         std::sin(2.0 * Pi * 6.0 * z / Lz + 0.12));
#endif

    // Sum-of-three triads has amplitude up to ~3; normalise back to ~1
    // so the final scaling is exactly `A` per component.
    const amrex::Real norm = 1.0 / 3.0;

    if (in_step_block) {
      AMREX_D_TERM(state(i, j, k, VELX) = 0.0;, state(i, j, k, VELY) = 0.0;
                   , state(i, j, k, VELZ) = 0.0;);
    } else {
      AMREX_D_TERM(state(i, j, k, VELX) = U + A * wy * pu * norm;
                   , state(i, j, k, VELY) = 0.0 + A * wy * pv * norm;
                   , state(i, j, k, VELZ) = 0.0 + A * wy * pw * norm;);
    }

    if (!(is_incompressible == 0)) {
      return;
    }

    state(i, j, k, TEMP) = prob_parm.T_mean;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
    state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

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
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    // Only the x_lo face is Inflow in this problem; idir=0, sgn=1.
    if (idir == 0 && sgn == 1) {
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      massfrac[O2_ID] = 0.233;
      massfrac[N2_ID] = 0.767;

      auto eos = pele::physics::PhysicsType::eos();

      s_ext[VELX] += prob_parm.meanFlowMag; // Add mean to normal flow only.

      s_ext[TEMP] = prob_parm.T_mean;

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
      }
    }
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal_eb(
    const amrex::Real /*xEBface*/[AMREX_SPACEDIM],
    const amrex::Real* /*xEBfaceNormal*/[AMREX_SPACEDIM],
    const amrex::CellData<const amrex::Real>& /*s_in*/,
    amrex::Real /*s_ext*/[NVAR],
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& /*prob_parm*/)
  {
    // Adiabatic wall — no explicit state override needed.
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bctype_eb(
    const amrex::Real /*xEBface*/[AMREX_SPACEDIM],
    int& EBflagType,
    amrex::Real& EBfaceFrac,
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& /*prob_parm*/)
  {
    EBflagType = pelelmex::BCTypeEB::wall_adiab;
    EBfaceFrac = 0.0; /* setting a safe value to ensure zero transport flux */
  }

  static void patchFlowVariables(
    const amrex::Geometry& /*geom*/,
    MyProbParm const& /*lprobparm*/,
    amrex::MultiFab& /*a_mf*/)
  {
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
