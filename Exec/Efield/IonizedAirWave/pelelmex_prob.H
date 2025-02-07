#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata(
  int i,
  int j,
  int k,
  int /*is_incompressible*/,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& /*aux*/,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

  AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
               , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
               , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real pert;
  amrex::Real y1, y2;

  pert = 0.0;
  if (prob_parm.pertmag > 0.0) {
    pert = prob_parm.pertmag *
           (1.0 * std::sin(2 * Pi * 4 * x / Lx) +
            1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) +
            0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) +
            1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) +
            0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx));
  }

  y1 = (y - prob_parm.standoff - 0.5 * dx[1] + pert) * 100;
  y2 = (y - prob_parm.standoff + 0.5 * dx[1] + pert) * 100;

  pele::physics::PMF::pmf(pmf_data, y1, y2, pmf_vals);

  state(i, j, k, TEMP) = pmf_vals[0];
  ;

  for (int n = 0; n < NUM_SPECIES; n++) {
    molefrac[n] = pmf_vals[3 + n];
  }
  eos.X2Y(molefrac, massfrac);

  AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
               , state(i, j, k, VELY) = pmf_vals[1] * 1e-2;
               , state(i, j, k, VELZ) = 0.0);

  amrex::Real P_cgs = prob_parm.P_mean * 10.0;

  // Density
  amrex::Real rho_cgs = 0.0;
  eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
  state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

  // Enthalpy
  amrex::Real h_cgs = 0.0;
  eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
  state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

  // Species mass
  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
  }

  amrex::Real yc = 0.5 * (prob_hi[1] + prob_lo[1]) + pert;
  state(i, j, k, NE) =
    2.0e13 * std::exp(-(y - yc) * (y - yc) / (Ly * Ly / 500.0));
  state(i, j, k, PHIV) = 0.0;
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
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == 1) {
    pele::physics::PMF::pmf(pmf_data, prob_lo[idir], prob_lo[idir], pmf_vals);

    AMREX_D_TERM(s_ext[VELX] = 0.0;, s_ext[VELY] = pmf_vals[1] * 1e-2;
                 , s_ext[VELZ] = 0.0);

    s_ext[TEMP] = pmf_vals[0];
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = pmf_vals[3 + n];
    }
    eos.X2Y(molefrac, massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
  }

  if (idir == 1) {  // y-dir
    if (sgn == 1) { // lo
      s_ext[NE] = 0.0;
      s_ext[PHIV] = prob_parm.phiV_loy;
    } else { // hi
      s_ext[NE] = 0.0;
      s_ext[PHIV] = prob_parm.phiV_hiy;
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& beta,
  amrex::GeometryData const& geomdata,
  amrex::Box const& domainBox,
  const int dir,
  const int beta_comp,
  const int nComp)
{
  amrex::ignore_unused(
    i, j, k, beta, geomdata, domainBox, dir, beta_comp, nComp);
  // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
  // otherwise this routine could be called for other face diffusivity (Temp,
  // velocity, ...)
}
#endif
