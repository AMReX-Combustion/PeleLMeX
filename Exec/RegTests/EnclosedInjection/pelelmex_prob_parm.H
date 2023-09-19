#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real T_mean = 300.0_rt;
  amrex::Real T_jet = 1000.0_rt;
  amrex::Real V_jet = 1.0_rt;
  amrex::Real jet_rad = 0.25_rt;
  amrex::Real bl_thickness = 0.01_rt;
};
#endif
