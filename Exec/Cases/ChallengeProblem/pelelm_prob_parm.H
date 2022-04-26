#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   // Chamber conditions
   amrex::Real P_mean        = 6079500.0;
   amrex::Real T_mean        = 900.0;
   amrex::Real Y_CH4_chamber = 0.020594;
   amrex::Real Y_O2_chamber  = 0.164305;

   // Injection parameters
   int nholes = 4;
   amrex::Real cone_angle = 45.0;
   amrex::Real centx      = 0.025;
   amrex::Real centy      = 0.025;
   amrex::Real r_circ     = 0.015;
   amrex::Real r_hole     = 0.00085;
   amrex::Real T_jet      = 470.0;
   amrex::Real vel_jet    = 28.0;
   amrex::Real inj_start  = 1e-6;
   amrex::Real inj_dur    = .5e-3;
   amrex::Real tau        = 0.01;
   amrex::Real Z          = 0.45;
};
#endif
