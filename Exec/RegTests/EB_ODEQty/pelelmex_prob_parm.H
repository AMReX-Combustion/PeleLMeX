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
  int meanFlowDir = 1;
  bool ode_qty_test = false;
  amrex::Real ode_IC = 0.0;
  amrex::Real ode_xy_lo = 0.0;
  amrex::Real ode_length = 0.0;
  amrex::Real ode_height = 0.0;
  amrex::Real ode_srcstrength = 0.0;
  bool composition_test = false;
  amrex::Real Y_CO2_0 = 0.0;
  amrex::Real Y_N2_0 = 0.25;
  amrex::Real Y_AR_0 = 0.75;
  amrex::Real extRhoYCO2 = 0.0;
  amrex::Real extRhoYCO2_ts = 0.0;
};
#endif
