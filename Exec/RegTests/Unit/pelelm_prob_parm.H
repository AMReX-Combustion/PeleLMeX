#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

 #include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real T_mean   = 300.0_rt;
};

#endif
