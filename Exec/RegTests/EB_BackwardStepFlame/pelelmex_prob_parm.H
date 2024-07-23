#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  // Shared params
  amrex::Real T_mean = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real Y_fuel = 0.0445_rt;
  amrex::Real Y_o2 = 0.223_rt;
  amrex::Real T_hot = 1800.0_rt;
  amrex::Real Twall = 300.0_rt;
  amrex::Real meanFlowMag = 0.0;

  bool ignite_flow = false;
  amrex::Real ignition_temperature = 1700.0;
};
#endif
