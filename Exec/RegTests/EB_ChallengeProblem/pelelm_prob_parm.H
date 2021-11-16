#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
    : amrex::Gpu::Managed
{
   amrex::Real P_mean   = 1013250.0_rt;
   amrex::Real T_mean   = 500.0_rt;
   amrex::Real rvort = 0.001;
   amrex::Real xvort = 0.0;
   amrex::Real yvort = 0.0;
   amrex::Real forcevort = 3.0;

   int nholes = 4;
   amrex::Real cone_angle = 20.0;
   amrex::Real centx = 0.025;
   amrex::Real centy = 0.025;
   amrex::Real r_circ = 0.03;
   amrex::Real r_hole = 0.004;
   amrex::Real T_jet = 300.0;
   amrex::Real vel_jet = 50.0;
};
#endif
