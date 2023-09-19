#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real inj_start = 1e-6_rt;
  amrex::Real inj_dur = .5e-3_rt;
  amrex::Real v_in = 28.0_rt;
  amrex::Real D = 1.7e-6_rt;
  amrex::Real Z = 0.45_rt;
  amrex::Real T_ox = 900.0_rt;
  amrex::Real T_fu = 470.0_rt;
  amrex::Real tau = 0.01_rt;

  amrex::Real Y_fuel[NUM_SPECIES] = {0.0};
  amrex::Real Y_ox[NUM_SPECIES] = {0.0};
  amrex::Real center_xy[2] = {0.0};
  amrex::Real H_fuel = {0.0};
  amrex::Real H_ox = {0.0};
};
#endif
