#ifndef PELELM_PROB_PARM_H
#define PELELM_PROB_PARM_H

#include <AMReX_REAL.H>

using namespace amrex::literals;

struct ProbParm
{
  amrex::Real P_mean = 101325.0_rt;
  amrex::Real standoff = 0.0_rt;
  amrex::Real pertmag = 0.0004_rt;
  amrex::Real Vin = -1.0;
  amrex::Real Vcoflow = 0.0;
  amrex::Real T_lips = 550.0;
  amrex::Real r_int = 0.0;
  amrex::Real r_ext = 0.0;
  amrex::Real splitx = 0.0;
  int phi = 1;
#ifdef PELE_USE_PLASMA
  amrex::Real electrode_radius = 0.0085;
  amrex::Real electrode_width = 0.0005;
  amrex::Real electrode_phiV = 6000.0;
  amrex::Real burner_phiV = 0.0;
  amrex::Real burner_Rext = 0.0;
#endif
};
#endif
