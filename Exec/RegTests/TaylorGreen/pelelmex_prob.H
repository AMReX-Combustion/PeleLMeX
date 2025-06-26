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
  amrex::Real T_mean = 300.0;
  amrex::Real reynolds = 1600.0;
  amrex::Real mach = 0.1;
  amrex::Real prandtl = 0.71;
  bool convecting = false;
  amrex::Real omega_x = 1.0;
  amrex::Real omega_y = 1.0;
  amrex::Real omega_z = 1.0;
  amrex::Real L = 0.0;
  amrex::Real rho0 = 0.0;
  amrex::Real v0 = 0.0;
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
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* dx = geomdata.CellSize();

    AMREX_D_TERM(
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
      , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    // TG functions
    amrex::Real u[3] = {0.0};
    u[0] = prob_parm.v0 * AMREX_D_TERM(
                            sin(prob_parm.omega_x * x / prob_parm.L),
                            *cos(prob_parm.omega_y * y / prob_parm.L),
                            *cos(prob_parm.omega_z * z / prob_parm.L));
    u[1] = -prob_parm.v0 * AMREX_D_TERM(
                             cos(prob_parm.omega_x * x / prob_parm.L),
                             *sin(prob_parm.omega_y * y / prob_parm.L),
                             *cos(prob_parm.omega_z * z / prob_parm.L));
    if (prob_parm.convecting) {
      u[0] += prob_parm.v0;
      u[1] += prob_parm.v0;
    }

    AMREX_D_TERM(state(i, j, k, VELX) = u[0];, state(i, j, k, VELY) = u[1];
                 , state(i, j, k, VELZ) = u[2]);

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
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
