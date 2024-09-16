#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <PelePhysics.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real standoff = 0.0_rt;
  amrex::Real pertmag = 0.0004_rt;
  amrex::Real pertlength = 0.008_rt;
  int meanFlowDir = 1;
  pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type> const*
    eosparm;
};
#endif
