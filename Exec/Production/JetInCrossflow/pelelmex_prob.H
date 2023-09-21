#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_Index.H>
#include <PelePhysics.H>
#include <pelelmex_prob_parm.H>
#include <PMFData.H>

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
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
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
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  molefrac[O2_ID] = prob_parm.X_O2;
  molefrac[N2_ID] = prob_parm.X_N2;
  eos.X2Y(molefrac, massfrac);

  state(i, j, k, TEMP) = prob_parm.ox_temp;

  // Set the initial velocity field
  amrex::Real x_velocity = 0.0;
  amrex::Real y_velocity = 0.0;
  amrex::Real z_velocity = 0.0;
  amrex::Real V_in = prob_parm.V_in;

#if (AMREX_SPACEDIM == 2)
  x_velocity = 0.0;
  y_velocity = 0.0;
#elif (AMREX_SPACEDIM == 3)
  x_velocity = V_in + prob_parm.pertmag_cf * V_in *
                        (sin(2 * Pi * 4 * x / Lx) * sin(2 * Pi * 5 * y / Ly) *
                           sin(2 * Pi * 3 * z / Lz) +
                         1.017 * sin(2 * Pi * 3 * x / Lx) *
                           sin(2 * Pi * 7 * y / Ly) * sin(2 * Pi * 5 * z / Lz));
  y_velocity = prob_parm.pertmag_cf * V_in *
               (sin(2 * Pi * 5 * x / Lx) * sin(2 * Pi * 3 * y / Ly) *
                  sin(2 * Pi * 7 * z / Lz) +
                1.017 * sin(2 * Pi * 9 * x / Lx) * sin(2 * Pi * 4 * y / Ly) *
                  sin(2 * Pi * 2 * z / Lz));
  z_velocity = prob_parm.pertmag_cf * V_in *
               (sin(2 * Pi * 2 * x / Lx) * sin(2 * Pi * 6 * y / Ly) *
                  sin(2 * Pi * 5 * z / Lz) +
                1.017 * sin(2 * Pi * 7 * x / Lx) * sin(2 * Pi * 4 * y / Ly) *
                  sin(2 * Pi * 3 * z / Lz));
#endif

  AMREX_D_TERM(state(i, j, k, VELX) = x_velocity;
               , state(i, j, k, VELY) = y_velocity;
               , state(i, j, k, VELZ) = z_velocity);

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
  const amrex::Real* dx = geomdata.CellSize();
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  amrex::Real V_in = prob_parm.V_in;
  amrex::Real jet_rad = prob_parm.jet_rad;
  int cf_dir = prob_parm.cf_dir;
  int jet_dir = prob_parm.jet_dir;

  amrex::Real X_air[NUM_SPECIES] = {0.0};
  amrex::Real X_jet[NUM_SPECIES] = {0.0};
  amrex::Real Y_air[NUM_SPECIES] = {0.0};
  amrex::Real Y_jet[NUM_SPECIES] = {0.0};
  amrex::Real U_air[AMREX_SPACEDIM] = {0.0};
  amrex::Real U_jet[AMREX_SPACEDIM] = {0.0};
  amrex::Real rho_air = 0.0;
  amrex::Real rho_jet = 0.0;
  amrex::Real radsq = 0.0;
  amrex::Real factor = 0.0;
  amrex::Real area = 0.0;
  amrex::Real jet_vel = 0.0;

  auto eos = pele::physics::PhysicsType::eos();

  // Find the directions
  int other_dir = 0;
  bool found_dir = false;

  while (!found_dir) {
    if (other_dir == jet_dir or other_dir == cf_dir) {
      other_dir += 1;
    } else {
      found_dir = true;
    }
  }

  // Initialize the jet and air molefracs
  X_air[O2_ID] = prob_parm.X_O2;
  X_air[N2_ID] = prob_parm.X_N2;
  X_jet[H2_ID] = prob_parm.jet_purity;
  X_jet[N2_ID] = 1.0 - prob_parm.jet_purity;

  // Get air/jet mass fraction and densities
  // Mass fraction
  eos.X2Y(X_air, Y_air);
  eos.X2Y(X_jet, Y_jet);
  // Density
  eos.PYT2R(prob_parm.P_mean * 10.0, Y_air, prob_parm.ox_temp, rho_air);
  eos.PYT2R(prob_parm.P_mean * 10.0, Y_jet, prob_parm.jet_temp, rho_jet);
  rho_air *= 1.0e3;
  rho_jet *= 1.0e3;

  // Initialize jet velocity BCs
#if (AMREX_SPACEDIM == 2)
  amrex::Abort("JIC not implemented in 2D");
#elif (AMREX_SPACEDIM == 3)
  area = (prob_hi[other_dir] - prob_lo[other_dir]) *
         (prob_hi[jet_dir] - prob_lo[jet_dir]);

  jet_vel = (prob_parm.global_eq_ratio * rho_air * area * V_in * Y_air[O2_ID]) /
            (8 * Pi * jet_rad * jet_rad * rho_jet * Y_jet[H2_ID]);

  if (cf_dir == 0) {
    U_air[0] = V_in;
  } else if (cf_dir == 1) {
    U_air[1] = V_in;
  } else {
    U_air[2] = V_in;
  }

  if (jet_dir == 0) {
    U_jet[0] = jet_vel;
  } else if (jet_dir == 1) {
    U_jet[1] = jet_vel;
  } else {
    U_jet[2] = jet_vel;
  }
#endif

  if (idir == jet_dir and (prob_parm.double_jet == 1 or sgn == 1)) {

    radsq = x[other_dir] * x[other_dir] + x[cf_dir] * x[cf_dir];
    factor =
      0.5 *
      (1.0 - tanh((sqrt(radsq) - jet_rad) / (prob_parm.bl_thickness / 4.0)));

    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = factor * X_jet[n] + (1.0 - factor) * X_air[n];
    }
    eos.X2Y(molefrac, massfrac);

    amrex::Real rho_cgs, P_cgs, RhoH_temp;
    P_cgs = prob_parm.P_mean * 10.0;

    s_ext[TEMP] =
      factor * prob_parm.jet_temp + (1.0 - factor) * prob_parm.ox_temp;

    eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
    s_ext[DENSITY] = rho_cgs * 1.0e3;

    eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY]; // CGS -> MKS conversion

    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = massfrac[n] * s_ext[DENSITY];
    }

    // Velocity
    amrex::Real time_factor = 0.0;
    if (time < prob_parm.init_time) {
      time_factor = (time / prob_parm.init_time) * (time / prob_parm.init_time);
    } else {
      time_factor = 1.0;
    }
    s_ext[VELX] =
      (1 - 0.5 * prob_parm.double_jet) * sgn * time_factor * factor * U_jet[0];
    s_ext[VELY] =
      (1 - 0.5 * prob_parm.double_jet) * sgn * time_factor * factor * U_jet[1];
    s_ext[VELZ] =
      (1 - 0.5 * prob_parm.double_jet) * sgn * time_factor * factor * U_jet[2];
  }

  if ((idir == cf_dir) and (sgn == 1)) {
    // Density
    s_ext[DENSITY] = rho_air;
    // Mass fraction
    for (int n = 0; n < NUM_SPECIES; n++) {
      s_ext[FIRSTSPEC + n] = Y_air[n] * s_ext[DENSITY];
    }
    // Temperature
    s_ext[TEMP] = prob_parm.ox_temp;
    amrex::Real RhoH_temp;
    // Enthalpy
    eos.TY2H(s_ext[TEMP], Y_air, RhoH_temp);
    s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];
    // Velocity
    s_ext[VELX] = U_air[0];
    s_ext[VELY] = U_air[1];
    s_ext[VELZ] = U_air[2];
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
