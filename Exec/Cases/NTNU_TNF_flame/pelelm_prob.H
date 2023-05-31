#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>



AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelm_initdata (int i, int j, int k,
                 int is_incompressible,
                 amrex::Array4<amrex::Real> const& state,
                 amrex::Array4<amrex::Real> const& aux,
                 amrex::GeometryData const& geomdata,
                 ProbParm const& prob_parm,
                 pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

    AMREX_D_TERM(,
                 const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1]);

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real pert = 0.0;
    amrex::Real y1, y2, z1, z2;

    amrex::Real massfrac[NUM_SPECIES]      = {0.0};
    amrex::Real molefrac[NUM_SPECIES]      = {0.0};
      
    state(i,j,k,TEMP) = prob_parm.T_in;

    for (int n = 0; n < NUM_SPECIES; n++)
      massfrac[n] = prob_parm.Ys[n];

    state(i,j,k,VELX) = 0;
    state(i,j,k,VELY) = 0;
    state(i,j,k,VELZ) = prob_parm.U_b;

    amrex::Real rho_cgs, P_cgs;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;            // CGS -> MKS conversion

    eos.TY2H(state(i,j,k,TEMP), massfrac, state(i,j,k,RHOH));
    state(i,j,k,RHOH) *= 1.0e-4 * state(i,j,k,DENSITY);   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
        state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();


  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  //load premixed composition to massfrac array
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = prob_parm.Ys[n];
  }


  auto eos = pele::physics::PhysicsType::eos();

  if(idir == 2 and sgn == 1){

    // amrex::GpuArray<amrex::Real, NUM_SPECIES+4> vel_from_file = {0.0};

    //radius = 0 corresponds to air passage inner wall
    amrex::Real r = sqrt(pow((0.0-x[0]),2)+pow((0.0-x[1]),2))-(prob_parm.bluff_inner_diam/2.0);

    //Compute y+ based on radius. Only compute it in the air passage region
    amrex::Real y_plus = 0.0;
    if(r >= 0.0 and r <= prob_parm.channel_height){
      //We have two boundary layers developing in the inner and outer walls
      if(r < prob_parm.channel_height/2.0){
        y_plus = prob_parm.friction_vel*r/prob_parm.mu;
      }
      else{
        y_plus = prob_parm.friction_vel*(prob_parm.channel_height - r)/prob_parm.mu;        
      }
    }

    //Compute u+: z-velocity profile
    //Transition from viscous sublayer to log layer is assumed to occur at y+=11
    amrex::Real u_plus = 0.0;
    if(y_plus <= 11.){
      u_plus = y_plus;
    }
    else{
      u_plus = prob_parm.one_over_vonKarman*log(y_plus) + prob_parm.B;
    }

    amrex::Real z_vel = prob_parm.friction_vel*u_plus;

    // s_ext[VELX] = 0.0;
    // s_ext[VELY] = 0.0;
    s_ext[VELZ] += z_vel;

    s_ext[TEMP] = prob_parm.T_in;

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    // printf("%e,%f \n",r,z_vel);

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
    }
  }
  else{

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    s_ext[VELX] = 0.0;
    s_ext[VELY] = 0.0;
    s_ext[VELZ] = 0.0;

    s_ext[TEMP] = prob_parm.T_in;

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int /*i*/, int /*j*/, int /*k*/,
           amrex::Array4<amrex::Real> const& /*beta*/,
           amrex::GeometryData const& /*geomdata*/,
           amrex::Box const& /*dBox*/,
           const int  /*dir*/,
           const int  /*state_comp*/,
           const int  /*nComp*/)
{
}
#endif
