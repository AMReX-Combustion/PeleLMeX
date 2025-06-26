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
#include <Utilities.H>

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

struct MyProbParm : public ProbParmDefault
{
  amrex::Real T0_gas = 500.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
  amrex::Real vel_gas = 0.;
  amrex::Real part_Re = -1.;
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
    state(i, j, k, TEMP) = prob_parm.T0_gas;
    state(i, j, k, VELX) = prob_parm.vel_gas;
#if (AMREX_SPACEDIM == 2)
    state(i, j, k, VELY) = 0.0;
#elif (AMREX_SPACEDIM == 3)
    state(i, j, k, VELY) = 0.0;
    state(i, j, k, VELZ) = 0.0;
#endif

    amrex::Real P_cgs = m2c::P(prob_parm.P_mean);

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = c2m::Rho(rho_cgs);

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
    state(i, j, k, RHOH) = c2m::H(h_cgs) * state(i, j, k, DENSITY);

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
      s_ext[VELX] = prob_parm.vel_gas;, s_ext[VELY] = 0.;, s_ext[VELZ] = 0.;);
    amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
    massfrac[N2_ID] = prob_parm.Y_N2;
    massfrac[O2_ID] = prob_parm.Y_O2;
    s_ext[TEMP] = prob_parm.T0_gas;
    amrex::Real rho_cgs, P_cgs;
    P_cgs = m2c::P(prob_parm.P_mean);

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, massfrac.data(), prob_parm.T0_gas, rho_cgs);

    s_ext[DENSITY] = c2m::Rho(rho_cgs);
    amrex::Real h_cgs;
    eos.TY2H(prob_parm.T0_gas, massfrac.data(), h_cgs);
    s_ext[RHOH] = c2m::H(h_cgs) * s_ext[DENSITY];
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
