#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  // Shared params
  int probType = 0;
  amrex::Real T_mean = 298.0_rt;
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real meanFlowMag = 0.0;
  int meanFlowDir = 1;

  // CoVo params
  amrex::Real rvort = 0.07;
  amrex::Real xvort = 0.5;
  amrex::Real yvort = 0.5;
  amrex::Real forcevort = 6.0;

  // CoGau & DiffGau params
  amrex::Real rgauss = 0.1;
  amrex::Real xgauss = 0.5;
  amrex::Real ygauss = 0.5;
  amrex::Real ampgauss = 0.1;
  int gauss_type = 0;

  amrex::Real gaussTime = 0.001;
  amrex::Real gaussDiff = 0.000001;
};
#endif
