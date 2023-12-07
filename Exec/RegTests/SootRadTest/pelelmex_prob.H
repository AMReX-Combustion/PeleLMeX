#ifndef PELELMEX_PROB_H
#define PELELMEX_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParmParse.H>

#include <pelelmex_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>

#include <PeleLMeX_Index.H>
#include <PelePhysics.H>
#include "SootModel.H"

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata(
  int i,
  int j,
  int k,
  int /*is_incompressible*/,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& aux,
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
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {{0.0}};
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real x1 = (x - prob_parm.standoff - 0.5 * dx[0]) * 100.;
  amrex::Real x2 = (x - prob_parm.standoff + 0.5 * dx[0]) * 100.;
  pele::physics::PMF::pmf(pmf_data, x1, x2, pmf_vals);
  state(i, j, k, TEMP) = pmf_vals[1];
  amrex::Real norm = 0.;
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = amrex::max(0., amrex::min(1., pmf_vals[3 + n]));
    norm += massfrac[n];
  }
  for (int n = 0; n < NUM_SPECIES; ++n) {
    massfrac[n] = massfrac[n] / norm;
  }
  AMREX_D_TERM(state(i, j, k, VELX) = pmf_vals[0] * 1.E-2;
               , state(i, j, k, VELY) = 0.;, state(i, j, k, VELZ) = 0.;);
  amrex::Real rho_cgs;
  auto P_cgs = prob_parm.P_mean * 10.;

  auto eos = pele::physics::PhysicsType::eos();
  eos.PYT2R(P_cgs, massfrac.data(), state(i, j, k, TEMP), rho_cgs);
  state(i, j, k, DENSITY) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(state(i, j, k, TEMP), massfrac.data(), state(i, j, k, RHOH));
  state(i, j, k, RHOH) = state(i, j, k, RHOH) * 1.0e-4 *
                         state(i, j, k, DENSITY); // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
  }
  for (int is = 0; is < NUM_SOOT_MOMENTS + 1; ++is) {
    state(i, j, k, FIRSTSOOT + is) = prob_parm.soot_vals[is];
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real x[AMREX_SPACEDIM],
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real time,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  if (sgn == 1) {
    pele::physics::PMF::pmf(pmf_data, prob_lo[idir], prob_lo[idir], pmf_vals);
    AMREX_D_TERM(s_ext[VELX] = pmf_vals[0] * 1.E-2;, s_ext[VELY] = 0.0;
                 , s_ext[VELZ] = 0.0;);

    s_ext[TEMP] = pmf_vals[1];

    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = pmf_vals[3 + n];
    }

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    auto eos = pele::physics::PhysicsType::eos();
    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }
    for (int is = 0; is < NUM_SOOT_MOMENTS + 1; ++is) {
      s_ext[FIRSTSOOT + is] = prob_parm.soot_vals[is];
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
