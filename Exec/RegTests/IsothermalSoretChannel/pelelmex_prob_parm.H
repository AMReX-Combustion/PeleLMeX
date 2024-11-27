#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <PeleLMeX_FlowControllerData.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real Vin = 10.0;
  amrex::Real Thigh = 700.0;
  amrex::Real Tlow = 300.0;
  amrex::Real phi = 1.0;

  int bathID{-1};
  int fuelID{-1};
  int oxidID{-1};
};
#endif
