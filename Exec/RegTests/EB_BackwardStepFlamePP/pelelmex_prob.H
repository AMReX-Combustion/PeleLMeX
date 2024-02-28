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
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* dx = geomdata.CellSize();

  const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
  const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];

  auto eos = pele::physics::PhysicsType::eos();

  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  massfrac[CH4_ID] = prob_parm.Y_fuel;
  massfrac[O2_ID] = prob_parm.Y_o2;
  massfrac[N2_ID] = 1.0 - massfrac[CH4_ID] - massfrac[O2_ID];

  AMREX_D_TERM(state(i, j, k, VELX) = prob_parm.meanFlowMag;
               , state(i, j, k, VELY) = 0.0;, state(i, j, k, VELZ) = 0.0);

  if (!(is_incompressible == 0)) {
    return;
  }

  if (x > 0.0 && y < 0.0) {
    state(i, j, k, TEMP) = prob_parm.T_hot;
    massfrac[CH4_ID] = 0.0;
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 1 - massfrac[CH4_ID] - massfrac[O2_ID];
  } else {
    state(i, j, k, TEMP) = prob_parm.T_mean;
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
  const int /*m_nAux*/,
  amrex::Real s_ext[NVAR],
  const int idir,
  const int sgn,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& prob_parm,
  pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
{
  if (sgn == 1) {
    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    auto eos = pele::physics::PhysicsType::eos();

    AMREX_D_TERM(s_ext[VELX] = prob_parm.meanFlowMag;, s_ext[VELY] = 0.0;
                 , s_ext[VELZ] = 0.0);

    massfrac[CH4_ID] = prob_parm.Y_fuel;
    massfrac[O2_ID] = prob_parm.Y_o2;
    massfrac[N2_ID] = 1.0 - massfrac[CH4_ID] - massfrac[O2_ID];

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

  // Set domain boundary Temp on the y-low face when using NoSlipWallIsotherm
  if (idir == 1 && sgn == 1) {
    s_ext[TEMP] = prob_parm.Twall;
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

struct EBhandle {
    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bctype_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
                   int & EBflagType,
                   amrex::GeometryData const& /*geomdata*/,
                   ProbParm const& /*prob_parm*/) const
    {
        if ( xEBface[0] > 0.0098 ) {
          EBflagType = pelelmex::BCTypeEB::wall_isot;
        } else {
          if (std::abs(xEBface[0]) <= 0.0015) {
            EBflagType = pelelmex::BCTypeEB::inflow;
          }
        }
    }

    AMREX_GPU_DEVICE
    AMREX_FORCE_INLINE
    void bcnormal_eb(const amrex::Real xEBface[AMREX_SPACEDIM],
                     const amrex::Real EBfaceNormal[AMREX_SPACEDIM], // Outward normal
                     const amrex::Real s_in[NVAR],
                     amrex::Real s_ext[NVAR],
                     const amrex::Real /*time*/,
                     amrex::GeometryData const& /*geomdata*/,
                     ProbParm const& prob_parm) const
    {
        // Isothermal wall
        s_ext[TEMP] = prob_parm.Twall;

        // Inflow
        if (std::abs(xEBface[0]) <= 0.0015) {
            auto eos = pele::physics::PhysicsType::eos();
            s_ext[VELX] = 0.0;
            s_ext[VELY] = prob_parm.EBinflow_vel;
            amrex::Real massfrac[NUM_SPECIES] = {0.0};
            massfrac[CH4_ID] = prob_parm.EBinflow_Yfuel;
            massfrac[O2_ID] = (1.0 - massfrac[CH4_ID]) * 0.233;
            massfrac[N2_ID] =  1.0 - massfrac[CH4_ID] - massfrac[O2_ID];

            amrex::Real P_cgs = prob_parm.P_mean * 10.0;

            s_ext[TEMP] = prob_parm.EBinflow_T;
            amrex::Real rho_cgs, RhoH_temp;
            eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
            s_ext[DENSITY] = rho_cgs * 1.0e3;

            eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
            s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

            for (int n = 0; n < NUM_SPECIES; n++) {
                s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
            }
        }
    }
};
#endif
