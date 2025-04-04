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
  pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type> const*
    eosparm;
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata(
    int i,
    int j,
    int k,
    int /*is_incompressible*/,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& aux,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();

    auto eos = pele::physics::PhysicsType::eos(prob_parm.eosparm);

    state(i, j, k, TEMP) = 300;
    amrex::Real P_cgs = prob_parm.P_mean * 10.0;
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    molefrac[N2_ID] = 0.79;
    molefrac[O2_ID] = 0.21;

    eos.X2Y(molefrac, massfrac);
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, FIRSTSPEC + n) = state(i, j, k, DENSITY) * massfrac[n];
    }
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
    state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
                 , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
                 , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    // want a cubic domain for ease
#if AMREX_SPACEDIM > 1
    AMREX_ALWAYS_ASSERT(
      AMREX_D_TERM(, amrex::almostEqual(Lx, Ly), &&amrex::almostEqual(Ly, Lz)));
#endif
    // find the centre of the domain
    const amrex::Real cen[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + 0.5 * Lx, prob_lo[1] + 0.5 * Ly, prob_lo[2] + 0.5 * Lz)};
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + (i + 0.5) * dx[0], prob_lo[1] + (j + 0.5) * dx[1],
      prob_lo[2] + (k + 0.5) * dx[2])};
    // radius from centre
    const amrex::Real r = std::sqrt(AMREX_D_TERM(
      (x[0] - cen[0]) * (x[0] - cen[0]), +(x[1] - cen[1]) * (x[1] - cen[1]),
      +(x[2] - cen[2]) * (x[2] - cen[2])));

    const amrex::Real sigma = Lx / 3.0;

    // impose radial Gaussian profile in all auxiliaries
    for (int n = 0; n < aux.nComp(); n++) {
      aux(i, j, k, n) = std::exp(-r * r / (sigma * sigma));
    }

    // constant upwards flow
    state(i, j, k, VELX) = 0.0;
#if AMREX_SPACEDIM == 2
    state(i, j, k, VELY) = 1.0;
#else
    state(i, j, k, VELY) = 0.0;
    state(i, j, k, VELZ) = 1.0;
#endif
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;

#endif
