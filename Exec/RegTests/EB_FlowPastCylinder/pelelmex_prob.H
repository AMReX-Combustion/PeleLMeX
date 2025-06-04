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
  amrex::Real meanFlowMag = 0.0;
  int meanFlowDir = 1;
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
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

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

    switch (prob_parm.meanFlowDir) {
    case 1:
      AMREX_D_TERM(s_ext[VELX] = prob_parm.meanFlowMag;, s_ext[VELY] = 0.0;
                   , s_ext[VELZ] = 0.0);
      break;
    case 2:
      AMREX_D_TERM(s_ext[VELY] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
                   , s_ext[VELZ] = 0.0);
      break;
    case 3:
      AMREX_D_TERM(s_ext[VELZ] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
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
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
