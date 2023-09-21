#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <pelelmex_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
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
  amrex::Array4<amrex::Real> const& aux,
  amrex::GeometryData const& geomdata,
  ProbParm const& prob_parm,
  PmfData const* /*pmf_data*/)
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

  AMREX_D_TERM(const amrex::Real x_c = prob_lo[0] + 0.5 * Lx;
               , const amrex::Real y_c = prob_lo[1] + 0.5 * Ly;
               , const amrex::Real z_c = prob_lo[2] + 0.5 * Lz);

  auto eos = pele::physics::PhysicsType::eos();
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = 0.233;
  massfrac[N2_ID] = 0.767;

  if (prob_parm.probType == 0) { // CoVo
    const amrex::Real deltax = x - prob_parm.xvort;
    const amrex::Real deltay = y - prob_parm.yvort;
    const amrex::Real d_sq = deltax * deltax + deltay * deltay;
    const amrex::Real r_sq = prob_parm.rvort * prob_parm.rvort;
    const amrex::Real u_vort =
      -prob_parm.forcevort * deltay / r_sq * exp(-d_sq / r_sq / 2.);
    const amrex::Real v_vort =
      prob_parm.forcevort * deltax / r_sq * exp(-d_sq / r_sq / 2.);
    const amrex::Real w_vort = 0.;

    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag + u_vort;
                 , state(i, j, k, VELY) = v_vort;, state(i, j, k, VELZ) = 0.0);

    if (is_incompressible)
      return;

    state(i, j, k, TEMP) = prob_parm.T_mean;
  } else if (prob_parm.probType == 1) { // CoGau

    AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
                 , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);

    if (is_incompressible)
      return;

    state(i, j, k, TEMP) = prob_parm.T_mean;
    const amrex::Real deltax = x - prob_parm.xgauss;
    const amrex::Real deltay = y - prob_parm.ygauss;
    const amrex::Real d_sq = deltax * deltax + deltay * deltay;
    const amrex::Real r_sq = prob_parm.rgauss * prob_parm.rgauss;
    if (prob_parm.gauss_type == 0) { // Spec
      massfrac[O2_ID] +=
        massfrac[O2_ID] * prob_parm.ampgauss * std::exp(-d_sq / r_sq);
      massfrac[N2_ID] = 1.0 - massfrac[O2_ID];
    } else if (prob_parm.gauss_type == 1) { // Temp
      state(i, j, k, TEMP) =
        prob_parm.T_mean * (1.0 + prob_parm.ampgauss * std::exp(-d_sq / r_sq));
    } else {
      amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
    }
  } else if (prob_parm.probType == 2) { // DifGau

    AMREX_D_TERM(state(i, j, k, VELX) = 0.0;, state(i, j, k, VELY) = 0.0;
                 , state(i, j, k, VELZ) = 0.0);

    if (is_incompressible)
      return;

    state(i, j, k, TEMP) = prob_parm.T_mean;
    const amrex::Real deltax = x - prob_parm.xgauss;
    const amrex::Real deltay = y - prob_parm.ygauss;
    const amrex::Real d_sq = deltax * deltax + deltay * deltay;
    const amrex::Real r_sq = prob_parm.rgauss * prob_parm.rgauss;
    if (prob_parm.gauss_type == 0) { // Spec
      massfrac[O2_ID] +=
        massfrac[O2_ID] * prob_parm.ampgauss * std::exp(-d_sq / r_sq);
      massfrac[N2_ID] = 1.0 - massfrac[O2_ID];
    } else if (prob_parm.gauss_type == 1) { // Temp
      state(i, j, k, TEMP) =
        prob_parm.T_mean * (1.0 + prob_parm.ampgauss * std::exp(-d_sq / r_sq));
    } else {
      amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
    }
  }

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
  PmfData const* /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  // amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  auto eos = pele::physics::PhysicsType::eos();
  if (sgn == 1) {
    // PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);

    s_ext[VELX] = 0.0;
#if (AMREX_SPACEDIM > 1)
    s_ext[VELY] = 0.0;
#if (AMREX_SPACEDIM > 2)
    s_ext[VELZ] = 0.0;
#endif
#endif

    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

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
}
#endif
