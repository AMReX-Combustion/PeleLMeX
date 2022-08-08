#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean        = 101325.0_rt;
   amrex::Real T_mean        = 300.0_rt;
   amrex::Real T_bubble      = 600.0_rt;
   amrex::Real bubble_rad    = 0.005_rt;
   amrex::Real bubble_y0     = 0.01;
   int         bubble_is_mix = 0;
   int           is_sym      = 0;
};
#endif
