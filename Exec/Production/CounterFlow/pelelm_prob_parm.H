#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T_ox = 298.0_rt;
  amrex::Real T_fuel = 298.0_rt;
  amrex::Real mdot = 0.01_rt;
  amrex::Real pertmag = 0.0004_rt;
  amrex::Real inertVel = 0.01_rt;

  amrex::Real jetRadius = 0.004_rt;
  amrex::Real inertRadius = 0.005_rt;

  int do_ignit = 0;
  amrex::Real ignitSphereRad = 0.003_rt;
  amrex::Real ignitT = 1800.0_rt;
};
#endif
