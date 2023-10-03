#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T_mean = 298.0_rt;
  amrex::Real rvort = 0.001;
  amrex::Real xvort = 0.0;
  amrex::Real yvort = 0.0;
  amrex::Real forcevort = 3.0;
};
#endif
