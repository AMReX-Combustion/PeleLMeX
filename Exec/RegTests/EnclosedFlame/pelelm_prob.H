#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
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
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data)
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

    AMREX_D_TERM(const amrex::Real xc = prob_lo[0];,// + 0.5 * Lx;,
                 const amrex::Real yc = prob_lo[1];,// + 0.5 * Ly;,
                 const amrex::Real zc = prob_lo[2];);// + 0.5 * Lz;);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert = 0.0;

    amrex::Real radius = std::sqrt(AMREX_D_TERM(  (x-xc)*(x-xc),
                                                + (y-yc)*(y-yc),
                                                + (z-zc)*(z-zc)));

    if (prob_parm.pertmag > 0.0) 
    {

#if ( AMREX_SPACEDIM == 2 )
       amrex::Real angle = std::atan2((y-yc),(x-xc));
       pert = prob_parm.pertmag * 
             (0.7 * std::sin(2*angle) +
              //1.023 * std::sin(2 * (angle - 0.004598) / radius2) +
              //0.945 * std::sin(3 * (angle - 0.00712435) / radius2) +
              0.9 * std::sin(8 * angle ) +
              1.4 * std::sin(4 * angle  ));
    }

    amrex::Real y1 = 3.6 - (radius - prob_parm.standoff - 1.42*dx[1] + pert)*100;
    amrex::Real y2 = 3.6 - (radius - prob_parm.standoff + 1.42*dx[1] + pert)*100;

#elif ( AMREX_SPACEDIM == 3 )

       amrex::Real angle = std::atan2((y-yc),(x-xc));
       amrex::Real angle2 = std::atan2((z-zc),(y-yc));
       pert = prob_parm.pertmag * 
             (0.5 * std::sin(2*angle) +
              0.7 * std::sin(8 * angle ) +
              1.1 * std::sin(4 * angle  )) *
             (1.04 * std::sin(2*angle2) +
              0.3 * std::sin(8 * angle2 ) +
              0.85 * std::sin(16 * angle2  ));
    }

    amrex::Real y1 = 3.6 - (radius - prob_parm.standoff - 1.42*dx[1] + pert)*100;
    amrex::Real y2 = 3.6 - (radius - prob_parm.standoff + 1.42*dx[1] + pert)*100;
#endif

    pele::physics::PMF::pmf(pmf_data,y1, y1, pmf_vals);

    temp(i,j,k) = pmf_vals[0];;

    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);

    vel(i,j,k,0) = 0.0;
#if ( AMREX_SPACEDIM == 2 )
    vel(i,j,k,1) = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 )
    vel(i,j,k,1) = 0.0;
    vel(i,j,k,2) = pmf_vals[1]*1e-2;
#endif

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
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
}
#endif
