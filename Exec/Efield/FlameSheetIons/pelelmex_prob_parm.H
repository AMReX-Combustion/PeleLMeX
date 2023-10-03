#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real standoff = 0.0_rt;
  amrex::Real pertmag = 0.0004_rt;
  int meanFlowDir = 1;

  amrex::Real phiV_loy = 0.0_rt;
  amrex::Real phiV_hiy = 0.0_rt;
};
#endif
