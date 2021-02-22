#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <EOS.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& vel,
                     amrex::Array4<amrex::Real> const& rho,
                     amrex::Array4<amrex::Real> const& rhoY,
                     amrex::Array4<amrex::Real> const& rhoH,
                     amrex::Array4<amrex::Real> const& temp,
                     amrex::Array4<amrex::Real> const& aux,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     PmfData const *pmf_data)
{
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx      = geomdata.CellSize();

    AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
                 const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
                 const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];); 

    AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                 const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                 const amrex::Real Lz = prob_hi[2] - prob_lo[2]); 

    AMREX_D_TERM(const amrex::Real x_c = prob_lo[0] + 0.5*Lx;,
                 const amrex::Real y_c = prob_lo[1] + 0.5*Ly;,
                 const amrex::Real z_c = prob_lo[2] + 0.5*Lz); 

    amrex::Real rad = std::sqrt(AMREX_D_TERM(  (x-x_c) * (x-x_c),
                                             + (y-y_c) * (y-y_c),
                                             + (z-z_c) * (z-z_c)));

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    AMREX_D_TERM(vel(i,j,k,0) = 0.0;,
                 vel(i,j,k,1) = 0.0;,
                 vel(i,j,k,2) = 0.0);
    
    if (is_incompressible) return;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    amrex::Real Gwidth = (AMREX_D_TERM(Lx,+Ly,+Lz))/(AMREX_SPACEDIM)/10.0;
    amrex::Real Gampl = 10.0;

    temp(i,j,k) = prob_parm.T_mean;
    temp(i,j,k) += Gampl * std::exp(-rad*rad/(2*Gwidth*Gwidth));

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    EOS::PYT2R(P_cgs, massfrac, temp(i,j,k), rho_cgs);
    rho(i,j,k) = rho_cgs * 1.0e3;

    // Enthalpy
    EOS::TY2H(temp(i,j,k), massfrac, rhoH(i,j,k));
    rhoH(i,j,k) *= 1.0e-4 * rho(i,j,k);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {  
       rhoY(i,j,k,n) = massfrac[n] * rho(i,j,k);
    }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int m_nAux,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  //amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  if (sgn == 1) {
    //PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);

    s_ext[VELX] = 0.0;
#if ( AMREX_SPACEDIM > 1 )
    s_ext[VELY] = 2.0;
#if ( AMREX_SPACEDIM > 2)
    s_ext[VELZ] = 0.0;
#endif
#endif

    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    s_ext[TEMP] = prob_parm.T_mean;

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    EOS::PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    EOS::TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
    }
  }
}
#endif
