#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
    amrex::Real T_mean = 300.0_rt;
    amrex::Real P_mean = 101325.0_rt;
    amrex::Real reynolds = 1600.0;
    amrex::Real mach = 0.1;
    amrex::Real prandtl = 0.71;
    bool convecting = false;
    amrex::Real omega_x = 1.0;
    amrex::Real omega_y = 1.0;
    amrex::Real omega_z = 1.0;
    amrex::Real L = 0.0;
    amrex::Real rho0 = 0.0;
    amrex::Real v0 = 0.0;
};
#endif
