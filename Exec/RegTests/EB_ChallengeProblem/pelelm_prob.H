#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <PelePhysics.H>

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
                     PmfData const * pmf_data)
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

    AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + 0.5 * Lx;,
                 const amrex::Real yc = prob_lo[1] + 0.5 * Ly;,
                 const amrex::Real zc = prob_lo[2] + 0.5 * Lz;);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    temp(i,j,k) = prob_parm.T_mean;

    amrex::Real r_xy = std::sqrt(x*x+y*y);
    amrex::Real angle = std::atan2(y,x);
    amrex::Real angular_vel = r_xy * 400.0 * std::tanh(300*(0.012-r_xy));

    vel(i,j,k,0) = angular_vel * std::sin(angle);
    vel(i,j,k,1) = -angular_vel * std::cos(angle);
    vel(i,j,k,2) = 0.0;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, temp(i,j,k), rho_cgs);
    rho(i,j,k) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(temp(i,j,k), massfrac, h_cgs);
    rhoH(i,j,k) = h_cgs * 1.0e-4 * rho(i,j,k);

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
  ProbParm const& prob_parm,
  PmfData const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
}
#endif
