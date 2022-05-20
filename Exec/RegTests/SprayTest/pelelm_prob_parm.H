#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
   amrex::Real P_mean   = 101325.0_rt;
   amrex::Real T0 = 500.;
   amrex::Real Y_O2 = 0.233;
   amrex::Real Y_N2 = 0.767;
   amrex::Real vel = 0.;
   amrex::IntVect partNum = amrex::IntVect(AMREX_D_DECL(100, 100, 100));
   amrex::Real partTemp = 300.;
   amrex::Real partDia = 1.E-3;
   amrex::RealVect partVel = amrex::RealVect(AMREX_D_DECL(0., 0., 0.));
   int numRedist = -1;
};
#endif
