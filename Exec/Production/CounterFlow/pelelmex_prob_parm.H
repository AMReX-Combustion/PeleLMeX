#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  // Gas params
  amrex::Real P_mean = 101325.0_rt;
  // amrex::Real T_ox     = 298.0_rt;
  // amrex::Real T_fuel   = 298.0_rt;
  // amrex::Real T_inert   = 298.0_rt;
  amrex::Real T_ox = 450.0_rt;
  amrex::Real T_fuel = 450.0_rt;
  amrex::Real T_inert = 450.0_rt;
  amrex::Real Y_O2_ox = 0.233_rt;
  amrex::Real Y_N2_ox = 0.767_rt;
  amrex::Real Y_N2_fuel = 0.0_rt;
  amrex::Real Y_fuel = 0.0_rt;
  amrex::Real Y_H2O_ign = 0.0_rt;
  amrex::Real Y_CO2_ign = 0.0_rt;
  amrex::Real mdot_ox = 0.01_rt;
  amrex::Real mdot_fuel = 0.01_rt;
  amrex::Real pertmag = 0.0004_rt;
  amrex::Real inertVel_ox = 0.01_rt;
  amrex::Real inertVel_fuel = 0.01_rt;
  amrex::Real jetRadius = 0.004_rt;
  amrex::Real inertRadius = 0.005_rt;

  // Misc
  int do_ignit = 0;
  amrex::Real ignitSphereRad = 0.003_rt;
  amrex::Real ignitT = 1800.0_rt;
};
#endif
