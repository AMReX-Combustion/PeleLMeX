#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include "PeleLMeX_Index.H"
#include "pelelmex_prob_parm.H"
#include "PMFData.H"
#include "PelePhysics.H"

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

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
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
  const amrex::Real y_low = prob_lo[1];

  const amrex::Real Lx = prob_hi[0] - prob_lo[0];
  const amrex::Real Ly = prob_hi[1] - prob_lo[1];

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  // Define air
  amrex::Real air_Y[NUM_SPECIES] = {0.0};
  air_Y[prob_parm.bathID] = 0.717;
  air_Y[prob_parm.oxidID] = 0.183;
  air_Y[prob_parm.fuelID] = 0.1;

  const amrex::Real T_low = 300.0;
  const amrex::Real T_high = 700.0;

  state(i, j, k, TEMP) = (y - y_low) / (Ly) * (T_high - T_low) + T_low;

  for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = air_Y[n];
  }

  state(i, j, k, VELX) = 10.0;
  state(i, j, k, VELY) = 0.0;

  amrex::Real rho_cgs, P_cgs;
  P_cgs = prob_parm.P_mean * 10.0;

  eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
  state(i, j, k, DENSITY) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(state(i, j, k, TEMP), massfrac, state(i, j, k, RHOH));
  state(i, j, k, RHOH) *=
    1.0e-4 * state(i, j, k, DENSITY); // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  if (sgn == 1) {
    s_ext[TEMP] = 300;
  } else if (sgn == -1) {
    s_ext[TEMP] = 700;
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
