#ifndef EBUSERKERNELS_H_
#define EBUSERKERNELS_H_

using namespace amrex;

#ifdef AMREX_USE_EB
#include "PelePhysics.H"
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBState(const amrex::Real xEBface[AMREX_SPACEDIM],
           amrex::Real s_ext[NVAR],
           const amrex::Real /*time*/,
           amrex::GeometryData const& /*geomdata*/,
           ProbParm const& prob_parm)
{
 // Isothermal wall
  s_ext[TEMP] = prob_parm.T_wall;

  amrex::Real r = sqrt(pow((0.0-xEBface[0]),2)+pow((0.0-xEBface[1]),2));

  // Inflow
  if (r > prob_parm.coflow_start_diam) {
      auto eos = pele::physics::PhysicsType::eos();
      s_ext[VELX] = 0.0;
      s_ext[VELY] = 0.0;
      s_ext[VELZ] = prob_parm.U_coflow;
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      massfrac[O2_ID] =  0.233;
      massfrac[N2_ID] =  0.767;

      amrex::Real P_cgs = prob_parm.P_mean * 10.0;

      s_ext[TEMP] = prob_parm.T_in;
      amrex::Real rho_cgs, RhoH_temp;
      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
          s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
setEBType(const amrex::Real xEBface[AMREX_SPACEDIM],
          amrex::Real & EBflagType,
          amrex::GeometryData const& /*geomdata*/,
          ProbParm const& /*prob_parm*/)
{
  if ( xEBface[0] > 0.0098 ) {
    EBflagType = 1.0;
  } else {
    EBflagType = 0.0;
  }
}

#endif
#endif