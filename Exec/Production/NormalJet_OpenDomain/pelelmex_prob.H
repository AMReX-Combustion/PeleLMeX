#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_func(
  amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau)
{
  amrex::Real eta_time;
  eta_time = 0.5 * (std::tanh((time - start) / (duration * tau)) *
                      std::tanh((-time + start + duration) / (duration * tau)) +
                    1.0);
  eta_time = std::min(std::max(eta_time, 0.0), 1.0);
  return eta_time;
}

struct MyProbParm : public ProbParmDefault
{
  amrex::Real inj_start = 1e-6;
  amrex::Real inj_dur = .5e-3;
  amrex::Real v_in = 28.0;
  amrex::Real D = 1.7e-6;
  amrex::Real Z = 0.45;
  amrex::Real T_ox = 900.0;
  amrex::Real T_fu = 470.0;
  amrex::Real tau = 0.01;

  amrex::Real Y_fuel[NUM_SPECIES] = {0.0};
  amrex::Real Y_ox[NUM_SPECIES] = {0.0};
  amrex::Real center_xy[2] = {0.0};
  amrex::Real H_fuel = {0.0};
  amrex::Real H_ox = {0.0};
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
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    state(i, j, k, TEMP) = prob_parm.T_ox;

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = prob_parm.Y_ox[n];
    }

    state(i, j, k, VELX) = 0.0;
    state(i, j, k, VELY) = 0.0;
    state(i, j, k, VELZ) = 1.0;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3; // CGS -> MKS conversion

    eos.TY2H(state(i, j, k, TEMP), massfrac, state(i, j, k, RHOH));
    state(i, j, k, RHOH) *=
      1.0e-4 * state(i, j, k, DENSITY); // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
    }
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real x[AMREX_SPACEDIM],
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int idir,
    const int sgn,
    const amrex::Real time,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real radius = {0.0};
    amrex::Real enthalpy;
    int ierr;

    if (sgn == 1 && idir == 2) {

      // Jet function
      radius = std::sqrt(
        std::pow(x[0] - prob_parm.center_xy[0], 2) +
        std::pow(x[1] - prob_parm.center_xy[1], 2));
      amrex::Real eta = 0.5 * (1.0 - std::tanh(
                                       2.0 * (radius - prob_parm.D * 0.5) /
                                       (0.1 * prob_parm.D * 0.5)));
      eta *=
        time_func(time, prob_parm.inj_start, prob_parm.inj_dur, prob_parm.tau);

      if (radius > prob_parm.D * 0.5) {
        s_ext[VELX] = 0.0;
        s_ext[VELY] = 0.0;
        s_ext[VELZ] = 0.0;
      }
      s_ext[VELZ] += eta * prob_parm.v_in;

      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] =
          eta * prob_parm.Y_fuel[n] + (1.0 - eta) * prob_parm.Y_ox[n];
      }
      enthalpy = eta * prob_parm.H_fuel + (1.0 - eta) * prob_parm.H_ox;

      GET_T_GIVEN_HY(enthalpy, massfrac, s_ext[TEMP], ierr);

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
