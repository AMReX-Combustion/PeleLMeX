#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>
#include <AMReX_GpuMemory.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 2431800.0_rt;
  amrex::Real V_in = 1.1e2_rt;
  amrex::Real jet_rad = 2.5e-4_rt;
  amrex::Real jet_temp = 5.5e2_rt;
  amrex::Real global_eq_ratio = 3.5e-1_rt;
  amrex::Real ox_temp = 7.5e2_rt;
  amrex::Real X_O2 = 2.1e-1_rt;
  amrex::Real X_N2 = 7.9e-1_rt;
  amrex::Real pertmag_cf = 0.0_rt;
  amrex::Real pertmag_jet = 0.0_rt;
  amrex::Real jet_purity = 0.995_rt;
  amrex::Real bl_thickness = 1.0e-4_rt;

  int init_time = -1;
  int double_jet = 1;
  int jet_dir = 0;
  int cf_dir = 2;
};
#endif
