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
  int is_incompressible,
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

  AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
               , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
               , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real u[3] = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real y1, y2;
  amrex::Real splitx = prob_lo[0];
  amrex::Real splity = prob_lo[1];
  amrex::Real rad_xy =
    std::sqrt((x - splitx) * (x - splitx) + (y - splity) * (y - splity));

  massfrac[N2_ID] = 0.767;
  massfrac[O2_ID] = 0.233;

  state(i, j, k, TEMP) = 298.0;

  if ((rad_xy) < (std::sqrt((z + 0.0009) / 0.01)) / 100.0) {
    amrex::Real pl1_angle = std::atan(0.01 / prob_parm.r_int);
    y1 = -(rad_xy - prob_parm.r_int) * std::cos(pl1_angle) +
         z * std::sin(pl1_angle);
    y2 = (rad_xy - prob_parm.r_int) * std::sin(pl1_angle) +
         z * std::cos(pl1_angle);
    y1 += prob_parm.standoff;
    y2 += prob_parm.standoff;

    pele::physics::PMF::pmf(pmf_data, y2, y2, pmf_vals);
    state(i, j, k, TEMP) = pmf_vals[0];
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac[n] = pmf_vals[3 + n];
    }
    // eos.X2Y(molefrac, massfrac);
  }

  AMREX_D_TERM(state(i, j, k, VELX) = 0.0;, state(i, j, k, VELY) = 0.1;
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
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();

  AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
               , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
               , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

  amrex::Real splitx = prob_lo[0];
  amrex::Real splity = prob_lo[1];
  amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
  massfrac_air[N2_ID] = 0.767;
  massfrac_air[O2_ID] = 0.233;

  if (prob_parm.phi == 1) { // Stoich
    massfrac_jet[CH4_ID] = 0.055;
  } else if (prob_parm.phi == 2) { // lean
    massfrac_jet[CH4_ID] = 0.0463;
  } else { // rich
    massfrac_jet[CH4_ID] = 0.0679;
  }
  massfrac_jet[O2_ID] = 0.233 * (1.0 - massfrac_jet[CH4_ID]);
  massfrac_jet[N2_ID] = 1.0 - massfrac_jet[CH4_ID] - massfrac_jet[O2_ID];

  auto eos = pele::physics::PhysicsType::eos();

  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = 0.0;
  }
  amrex::Real rad_xy = std::sqrt(
    (x[0] - splitx) * (x[0] - splitx) + (x[1] - splity) * (x[1] - splity));

  // Get zone
  int zone = -1;
  if (rad_xy < prob_parm.r_int) {
    zone = 3; // jet
  } else {
    zone = 1; // coflow
  }

  if (idir == 2 && sgn == 1) {
    s_ext[VELX] = 0.0;
    s_ext[VELY] = 0.0;
    if (zone == 3) {
      s_ext[VELZ] =
        2.0 * prob_parm.Vin *
        (1.0 - (rad_xy / prob_parm.r_int) * (rad_xy / prob_parm.r_int));
    } else if (zone == 1) {
      s_ext[VELZ] =
        prob_parm.Vcoflow *
        (1.0 + 1.5 * (rad_xy - prob_parm.r_ext) / (0.025 - prob_parm.r_ext));
      if (rad_xy <= prob_parm.r_ext)
        s_ext[VELZ] = 0.0;
    }

    s_ext[TEMP] = 298.0;

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;
    if (zone == 1) {
      eos.PYT2R(P_cgs, massfrac_air, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac_air, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac_air[n] * s_ext[DENSITY];
      }
      if (rad_xy <= prob_parm.r_ext)
        s_ext[TEMP] = prob_parm.T_lips;
    } else if (zone == 3) {
      s_ext[TEMP] = 302.0; // Jet is slightly warmed up
      eos.PYT2R(P_cgs, massfrac_jet, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac_jet, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac_jet[n] * s_ext[DENSITY];
      }
    }
  }

#ifdef PELE_USE_PLASMA
  if (idir == 2 && sgn == 1) { // Burner side
    amrex::Real dampWidth = 0.0005;
    amrex::Real lipFunction =
      0.5 * (std::tanh((rad_xy - prob_parm.r_int) / dampWidth) -
             std::tanh((rad_xy - prob_parm.burner_Rext) / dampWidth));
    s_ext[PHIV] =
      0.01 * prob_parm.burner_phiV + lipFunction * prob_parm.burner_phiV;
  } else if (idir == 2 && sgn == -1) { // electrode side
    amrex::Real dampWidth = prob_parm.electrode_width;
    amrex::Real lipFunction =
      0.5 *
      (std::tanh((rad_xy - (prob_parm.electrode_radius - 0.0025)) / dampWidth) -
       std::tanh((rad_xy - (prob_parm.electrode_radius + 0.0025)) / dampWidth));
    s_ext[PHIV] =
      0.01 * prob_parm.electrode_phiV + lipFunction * prob_parm.electrode_phiV;
  }
#endif
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
  // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
  // otherwise this routine could be called for other face diffusivity (Temp,
  // velocity, ...)
  if (beta_comp == 0 && nComp == NUM_SPECIES) {
    if (dir == 2 && k <= domainBox.smallEnd(dir)) {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      const amrex::Real* prob_hi = geomdata.ProbHi();
      const amrex::Real* dx = geomdata.CellSize();

      AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
                   , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
                   , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);
      amrex::Real rad_xy = std::sqrt(x * x + y * y);
      if (rad_xy >= 0.0065 && rad_xy <= 0.0075) { // hard coded size of the
                                                  // burner lip ...
        for (int n = 0; n < nComp; n++) {
          beta(i, j, k, n) = 0.0;
        }
      }
    }
  }
}
#endif
