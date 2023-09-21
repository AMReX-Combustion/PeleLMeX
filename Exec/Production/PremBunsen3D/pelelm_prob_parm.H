#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real standoff = 0.0_rt;
  amrex::Real pertmag = 0.0004_rt;
  amrex::Real Vin = -1.0;
  amrex::Real Vcoflow = 0.0;
  amrex::Real slot_width = 0.0;
  amrex::Real splitx = 0.0;
};
#endif
