#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

 #include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real T_mean   = 300.0_rt;
   amrex::Real meanFlowMag = 0.0;
   int  meanFlowDir = 1;

   // CoVo params
   amrex::Real rvort = 0.07;
   amrex::Real xvort = 0.5;
   amrex::Real yvort = 0.5;
   amrex::Real forcevort = 6.0;
};

#endif
