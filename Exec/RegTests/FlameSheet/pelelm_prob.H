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

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert = 0.0;

    if (prob_parm.pertmag > 0.0) 
    {

#if ( AMREX_SPACEDIM == 2 )
       pert = prob_parm.pertmag * 
             (1.0 * std::sin(2 * Pi * 4 * x / Lx) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx));
    }

    amrex::Real y1 = (y - prob_parm.standoff - 0.5*dx[1] + pert)*100;
    amrex::Real y2 = (y - prob_parm.standoff + 0.5*dx[1] + pert)*100;

#elif ( AMREX_SPACEDIM == 3 )

       pert = prob_parm.pertmag *
             (1.0 * std::sin(2 * Pi * 4 * x / Lx) *
                std::sin(2 * Pi * 5 * y / Ly) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) *
                std::sin(2 * Pi * 4 * (y - 0.0053765) / Ly) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) *
                std::sin(2 * Pi * 3 * (y - 0.02137) / Ly) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) *
                std::sin(2 * Pi * 6 * (y - 0.018) / Ly) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx));
    }

    amrex::Real y1 = (z - prob_parm.standoff - 0.5*dx[2] + pert)*100;
    amrex::Real y2 = (z - prob_parm.standoff + 0.5*dx[2] + pert)*100;
#endif


    PMF::pmf(pmf_data,y1, y2, pmf_vals);

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
  PmfData const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == 1) {
    PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);

    s_ext[VELX] = 0.0;
#if ( AMREX_SPACEDIM == 2 )
    s_ext[VELY] = pmf_vals[1]*1e-2;
#elif ( AMREX_SPACEDIM == 3 )
    s_ext[VELY] = 0.0;
    s_ext[VELZ] = pmf_vals[1]*1e-2;
#endif

    s_ext[TEMP] = pmf_vals[0];
    for (int n = 0; n < NUM_SPECIES; n++){
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
    }
  }
}
#endif
