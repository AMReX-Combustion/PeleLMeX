#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

struct MyProbParm : public ProbParmDefault
{
  amrex::Real T0 = 500.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  amrex::Real vel = 0.;
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
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[N2_ID] = prob_parm.Y_N2;
    massfrac[O2_ID] = prob_parm.Y_O2;
    state(i, j, k, TEMP) = prob_parm.T0;

    state(i, j, k, VELX) = prob_parm.vel;
#if (AMREX_SPACEDIM == 2)
    state(i, j, k, VELY) = 0.0;
#elif (AMREX_SPACEDIM == 3)
    state(i, j, k, VELY) = 0.0;
    state(i, j, k, VELZ) = 0.0;
#endif

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
    AMREX_D_TERM(
      s_ext[VELX] = prob_parm.vel;, s_ext[VELY] = 0.;, s_ext[VELZ] = 0.;);
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
    massfrac[N2_ID] = prob_parm.Y_N2;
    massfrac[O2_ID] = prob_parm.Y_O2;
    s_ext[TEMP] = prob_parm.T0;
    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, massfrac.data(), prob_parm.T0, rho_cgs);
    amrex::Real rho_mks = rho_cgs * 1.0e3; // CGS -> MKS conversion
    s_ext[DENSITY] = rho_mks;
    amrex::Real rhoH;
    eos.TY2H(prob_parm.T0, massfrac.data(), rhoH);
    s_ext[RHOH] = rhoH * 1.0e-4 * rho_mks;
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * rho_mks;
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
