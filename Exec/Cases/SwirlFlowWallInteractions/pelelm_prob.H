#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
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
                 pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
    auto eos = pele::physics::PhysicsType::eos();

    amrex::Real massfrac_eq[NUM_SPECIES] = {0.0};

    amrex::Real y_lo = 0.0001;
    amrex::Real Gauss_maxT = 1800.0;
    amrex::Real Gauss_T_width = 0.0024;
    amrex::Real Gauss_Spec_width = 0.0004;


    amrex::Real T_sigma = 0.002;
    amrex::Real Spec_sigma = 0.0001;
    amrex::Real T_max_exponential = 2000;


    amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    amrex::Real Ro_swirler = 0.0162;
    amrex::Real Ri_swirler = 0.0101;
    amrex::Real R = Ri_swirler/Ro_swirler;
    amrex::Real x0 = 0.0;
    amrex::Real y0 = 0.0;
    amrex::Real H = 0.00001;
    amrex::Real Do_swirler = 0.0324;
    amrex::Real Di_swirler = 0.0202;
    amrex::Real D_coflow = Do_swirler;
    amrex::Real Ri_coflow = Ro_swirler;
    amrex::Real c1 = -0.0067;
    amrex::Real c2 = 3.39e-7;
    amrex::Real c3 = -21.081;
    amrex::Real A = 0.0096;
    amrex::Real B = -3.869e-11;
    amrex::Real F_theta = -7.1884e5;


    massfrac_air[prob_parm.oxidID] = 0.233;  
    massfrac_air[prob_parm.bathID] = 0.767;

    state(i,j,k,TEMP) = 300.0;

    for (int n = 0; n < NUM_SPECIES; n++)
    {
        massfrac[n] = massfrac_air[n];
    }

    state(i,j,k,VELX) = 0.0;
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = 1.0;

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

  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};  
  amrex::Real massfrac[NUM_SPECIES] = {0.0}; 

  amrex::Real Ro_swirler = 0.0162;
  amrex::Real Ri_swirler = 0.0101;
  amrex::Real R = Ri_swirler/Ro_swirler;
  amrex::Real D_jet = 0.02;
  amrex::Real x0 = 0.0;
  amrex::Real y0 = 0.0;
  amrex::Real H = 0.00001;
  amrex::Real Do_swirler = 0.0324;
  amrex::Real Di_swirler = 0.0202;
  amrex::Real D_coflow = Do_swirler;
  amrex::Real Ro_coflow = prob_hi[0] ;
  amrex::Real Ri_coflow = Ro_swirler;
  amrex::Real c1 = -0.0067;
  amrex::Real c2 = 3.39e-7;
  amrex::Real c3 = -21.081;
  amrex::Real A = 0.0096;
  amrex::Real B = -3.869e-11;
  amrex::Real F_theta = -7.1884e5;

  massfrac_air[prob_parm.oxidID] = 0.233;
  massfrac_air[prob_parm.bathID] = 0.767;

  amrex::Real vel_swirler_x = 3 * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*((x[1]/(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))*(F_theta)*(((std::pow((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))),2))/3)+(c1*(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))+((c2)/((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2)))))));

  amrex::Real vel_swirler_y = 3 * 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*(x[0]/(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))*(F_theta)*(((std::pow((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))),2))/3)+(c1*(std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))+((c2)/((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2))))));


  s_ext[VELX] = -vel_swirler_x;
  s_ext[VELY] = vel_swirler_y;

  amrex::Real vel_swirler_z = 0.25*(1-((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Do_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Do_swirler/2))/(H))))*(1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(Di_swirler/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(Di_swirler/2))/(H))))*prob_parm.U_b*(c3)*((std::pow(((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2)))/(Ro_swirler)),2))-(((std::pow(R,2)-1)/(std::log(R)))*std::log((std::sqrt(std::pow(x[0],2)+std::pow(x[1],2)))/Ro_swirler))-1);

  amrex::Real vel_coflow_z = prob_parm.U_c * 0.5 * (1+((std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_coflow/2))/(H))))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_coflow/2))/(H))));

  amrex::Real vel_jet = prob_parm.U_b*( 1 - (0.5*(1+(std::tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))-(D_jet/2))/(H)))*(tanh((std::sqrt(std::pow(x[0]-x0,2)+std::pow(x[1]-y0,2))+(D_jet/2))/(H))))));

  s_ext[VELZ] = vel_swirler_z + vel_coflow_z + vel_jet;

  s_ext[TEMP] = prob_parm.T_in;

  for (int n = 0; n < NUM_SPECIES; ++n) {
     massfrac[n] = massfrac_air[n];
  }

  amrex::Real rho_cgs, P_cgs, RhoH_temp;
  P_cgs = prob_parm.P_mean * 10.0;

  auto eos = pele::physics::PhysicsType::eos();
  
  eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
  s_ext[DENSITY] = rho_cgs * 1.0e3;

  eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
  s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
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
