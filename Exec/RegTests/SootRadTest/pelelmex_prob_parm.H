#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::GpuArray<amrex::Real, NUM_SOOT_MOMENTS + 1> soot_vals;
  amrex::Real P_mean = 101325.0;
  amrex::Real standoff = 0.;
};

#endif
