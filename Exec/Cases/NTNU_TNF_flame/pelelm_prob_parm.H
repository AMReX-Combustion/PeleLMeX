#ifndef PELELM_PROB_PARM_H_
#define PELELM_PROB_PARM_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>
#include <PelePhysics.H>

using namespace amrex::literals;

struct ProbParm
{

  amrex::Real P_mean = 101325.0;
  amrex::Real T_in   = 300.0;
  amrex::Real T_wall  = 300.0_rt;
  amrex::Real U_b    = 5.0;

  //Crossflow related parameters
  amrex::Real U_crossflow   = 0.0;
  amrex::Real z_height_tanh = 0.0;
  amrex::Real delta         = 1.0;

  amrex::Real U_coflow   = 0.0;
  amrex::Real coflow_start_diam = 1.0e+10;

  amrex::Real one_over_vonKarman = 0.0; // 1/kappa -> one over von Karman constant
  amrex::Real B = 0.0; //constant used to compute turbulent log layer
  amrex::Real bluff_inner_diam = 0.0;
  amrex::Real channel_height   = 0.0;
  amrex::Real phi_in           = 0.0;
  amrex::Real mu               = 0.0;
  amrex::Real friction_vel     = 0.0;
  amrex::Real moles_NH3        = 0.0;
  amrex::Real moles_H2         = 0.0;
  amrex::Real moles_N2         = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Ys = {{0.0}};

  bool ignition = false;
  amrex::GpuArray<amrex::Real, 3> ign_region;
  amrex::Real T0 = 900.;
  amrex::Real Y_O2 = 0.233;
  amrex::Real Y_N2 = 0.767;
};
#endif
