#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm : amrex::Gpu::Managed
{
  // Shared params
  amrex::Real T_mean = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real meanFlowMag = 0.0;
  amrex::Real ode_IC = 10.0;
  amrex::Real ode_xy_lo = -0.01;
  amrex::Real ode_length = 0.01;
  amrex::Real ode_height = 0.020;
  amrex::Real ode_srcstrength = -10.0;
  int meanFlowDir = 1;
};
#endif
