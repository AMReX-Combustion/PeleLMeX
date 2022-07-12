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

   // Simplification
   int doQuarterDomain    = 0;

   // HIT init
   int hitIC = 0;
   int input_resolution = 0;
   amrex::Real uin_norm = 1.0;
   amrex::Real lscale = 1.0;
   amrex::Real offset = 0.0;
   amrex::Real urms0 = 0.0;
   amrex::Real Linput = 0.0;
   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_lo = {0.0, 0.0, 0.0};
   amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> win_hi = {0.0, 0.0, 0.0};
   amrex::Real win_slope = 1.0;
   amrex::Real* d_uinput = nullptr;
   amrex::Real* d_vinput = nullptr;
   amrex::Real* d_winput = nullptr;
   amrex::Real* d_xarray = nullptr;
   amrex::Real* d_xdiff  = nullptr;
};
#endif
