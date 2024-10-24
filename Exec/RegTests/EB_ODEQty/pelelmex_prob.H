#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pelelmex_initdata(
  int i,
  int j,
  int k,
  int is_incompressible,
  amrex::Array4<amrex::Real> const& state,
  amrex::Array4<amrex::Real> const& /*aux*/,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = 0.233;
  massfrac[N2_ID] = 0.767;

  switch (prob_parm.meanFlowDir) {
  case 1:
    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
    break;
  case -1:
    AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);
    break;
  case 2:
    AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                 , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  case -2:
    AMREX_D_TERM(state(i, j, k, VELX) = 0.0;
                 , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  case 3:
    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  case -3:
    AMREX_D_TERM(state(i, j, k, VELX) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = -prob_parm.meanFlowMag;
                 , state(i, j, k, VELZ) = 0.0);
    break;
  }

  if (!(is_incompressible == 0)) {
    return;
  }

  state(i, j, k, TEMP) = prob_parm.T_mean;

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

  // Initial ODE quantity
#if NUM_ODE > 0
  amrex::Real ode_qty[NUM_ODE] = {0.0};

  // Initial state for ode quantities
  amrex::Real ode_qty_0[NUM_ODE] = {0.0};
  for (int n = 0; n < NUM_ODE; n++) {
    ode_qty_0[n] = (n + 1) * prob_parm.ode_IC;
  }

  // Current x,y,z locations
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();
  const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
    prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
    prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
    prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2])};

  // Create NUM_ODE boxes from (x_strt,y_strt) to (x_end,y_end)
  for (int n = 0; n < NUM_ODE; n++) {
    amrex::Real x_strt = prob_parm.ode_xy_lo + n * 4. * prob_parm.ode_length;
    amrex::Real x_end = x_strt + prob_parm.ode_length;
    amrex::Real y_strt = prob_parm.ode_xy_lo;
    amrex::Real y_end = y_strt + prob_parm.ode_height;
    if ((x[0] > x_strt) && (x[1] > y_strt)) {
      if ((x[0] < x_end) && (x[1] < y_end)) {
        ode_qty[n] = ode_qty_0[n];
      }
    }
  }
  // Add to state variable
  for (int n = 0; n < NUM_ODE; n++) {
    state(i, j, k, FIRSTODE + n) = ode_qty[n];
  }

#endif
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
bcnormal(
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int /*m_node*/,
  amrex::Real s_ext[NVAR],
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  auto eos = pele::physics::PhysicsType::eos();

  switch (prob_parm.meanFlowDir) {
  case 1:
    AMREX_D_TERM(s_ext[VELX] = prob_parm.meanFlowMag;, s_ext[VELY] = 0.0;
                 , s_ext[VELZ] = 0.0);
    break;
  case 2:
    AMREX_D_TERM(s_ext[VELY] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
                 , s_ext[VELZ] = 0.0);
    break;
  case 3:
    AMREX_D_TERM(s_ext[VELZ] = prob_parm.meanFlowMag;, s_ext[VELX] = 0.0;
                 , s_ext[VELY] = 0.0);
    break;
  }

  massfrac[O2_ID] = 0.333;
  massfrac[N2_ID] = 0.667;

  s_ext[TEMP] = prob_parm.T_mean;

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
