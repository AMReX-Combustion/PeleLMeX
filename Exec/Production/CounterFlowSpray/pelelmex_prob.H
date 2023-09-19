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

  AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + Lx / 2.0;
               , const amrex::Real yc = prob_lo[1] + Ly / 2.0;
               , const amrex::Real zc = prob_lo[2] + Lz / 2.0);

  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_cold[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_hot[NUM_SPECIES] = {0.0};

  massfrac_cold[NXC7H16_ID] = 0.05;
  massfrac_cold[O2_ID] = 0.233 * (1.0 - massfrac_cold[NXC7H16_ID]);
  ;
  massfrac_cold[N2_ID] = 1.0 - massfrac_cold[NXC7H16_ID] - massfrac_cold[O2_ID];

  state(i, j, k, TEMP) = 298.0;

  if (prob_parm.do_ignit) {
    amrex::Real radiusSq = AMREX_D_TERM(
      (x - xc) * (x - xc), +(y - yc) * (y - yc), +(z - zc) * (z - zc));
    amrex::Real radius = std::sqrt(radiusSq);
    amrex::Real mixingWidth = 0.1 * prob_parm.ignitSphereRad;
    amrex::Real mixingFunction =
      0.5 *
      (1.0 + std::tanh((prob_parm.ignitSphereRad - radius) / mixingWidth));
    state(i, j, k, TEMP) =
      mixingFunction * prob_parm.ignitT + (1.0 - mixingFunction) * 298.0;
    massfrac_hot[NXC7H16_ID] = 0.00;
    massfrac_hot[O2_ID] = 0.233 * (1.0 - massfrac[NXC7H16_ID]);
    massfrac_hot[N2_ID] = 1.0 - massfrac_hot[NXC7H16_ID] - massfrac_hot[O2_ID];

    massfrac[NXC7H16_ID] = mixingFunction * massfrac_hot[NXC7H16_ID] +
                           (1.0 - mixingFunction) * massfrac_cold[NXC7H16_ID];
    massfrac[O2_ID] = mixingFunction * massfrac_hot[O2_ID] +
                      (1.0 - mixingFunction) * massfrac_cold[O2_ID];
    massfrac[N2_ID] = mixingFunction * massfrac_hot[N2_ID] +
                      (1.0 - mixingFunction) * massfrac_cold[N2_ID];
  } else {
    massfrac[NXC7H16_ID] = massfrac_cold[NXC7H16_ID];
    massfrac[O2_ID] = massfrac_cold[O2_ID];
    massfrac[N2_ID] = massfrac_cold[N2_ID];
  }

  AMREX_D_TERM(state(i, j, k, VELX) = 0.0;, state(i, j, k, VELY) = 0.0;
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

  amrex::Real splitx = prob_lo[0] + 0.5 * Lx;
  amrex::Real splity = prob_lo[1] + 0.5 * Ly;

  amrex::Real massfrac_ox[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_inert[NUM_SPECIES] = {0.0};
  amrex::Real massfrac_mix[NUM_SPECIES] = {0.0};

  massfrac_ox[N2_ID] = 0.767;
  massfrac_ox[O2_ID] = 0.233;
  massfrac_fuel[NXC7H16_ID] = 0.0;
  massfrac_fuel[N2_ID] = 1.0;
  massfrac_inert[N2_ID] = 1.0;

  auto eos = pele::physics::PhysicsType::eos();

  for (int n = 0; n < NVAR; n++) {
    s_ext[n] = 0.0;
  }

  // Get zone
  int zone = -1;
  if (x[0] < splitx) {
    zone = 1; // Oxidizer side jet
  } else {
    zone = 2; // Fuel side jet
  }
  amrex::Real radius = std::abs(x[1] - splity);
  amrex::Real mixingWidth = 0.1 * prob_parm.jetRadius;
  amrex::Real mixingFunction =
    0.5 * (1.0 + std::tanh((prob_parm.jetRadius - radius) / mixingWidth));

  amrex::Real P_cgs = prob_parm.P_mean * 10.0;
  amrex::Real rho_cgs, RhoH_temp;
  amrex::Real rhoOx = -1.0;
  {
    eos.PYT2R(P_cgs, massfrac_ox, prob_parm.T_ox, rhoOx);
    rhoOx *= 1.0e3;
  }
  amrex::Real rhoFuel = -1.0;
  {
    eos.PYT2R(P_cgs, massfrac_fuel, prob_parm.T_fuel, rhoFuel);
    rhoFuel *= 1.0e3;
  }
  if (zone == 1) {
    s_ext[TEMP] =
      prob_parm.T_ox * mixingFunction + 298.0 * (1.0 - mixingFunction);
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac_mix[n] = massfrac_ox[n] * mixingFunction +
                        massfrac_inert[n] * (1.0 - mixingFunction);
    }
    eos.PYT2R(P_cgs, massfrac_mix, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac_mix, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac_mix[n] * s_ext[DENSITY];
    }
  } else if (zone == 2) {
    s_ext[TEMP] =
      prob_parm.T_fuel * mixingFunction + 298.0 * (1.0 - mixingFunction);
    for (int n = 0; n < NUM_SPECIES; n++) {
      massfrac_mix[n] = massfrac_fuel[n] * mixingFunction +
                        massfrac_inert[n] * (1.0 - mixingFunction);
    }
    eos.PYT2R(P_cgs, massfrac_mix, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac_mix, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac_mix[n] * s_ext[DENSITY];
    }
  }

  s_ext[VELY] = 0.0;
  if (zone == 1) {
    amrex::Real vel =
      (prob_parm.mdot / rhoOx) * mixingFunction +
      prob_parm.inertVel /* (1.0 + (radius/0.025))*/ * (1.0 - mixingFunction);
    s_ext[VELX] = vel;
  } else if (zone == 2) {
    amrex::Real vel =
      (prob_parm.mdot / rhoFuel) * mixingFunction +
      prob_parm.inertVel /* (1.0 + (radius/0.025))*/ * (1.0 - mixingFunction);
    s_ext[VELX] = -vel;
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
