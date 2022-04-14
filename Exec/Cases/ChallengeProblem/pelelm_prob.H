#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void pelelm_initdata(int i, int j, int k,
                     int is_incompressible,
                     amrex::Array4<amrex::Real> const& state,
                     amrex::Array4<amrex::Real> const& /*aux*/,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     pele::physics::PMF::PmfData::DataContainer const * pmf_data)
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

    AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + 0.5 * Lx;,
                 const amrex::Real yc = prob_lo[1] + 0.5 * Ly;,
                 const amrex::Real zc = prob_lo[2] + 0.5 * Lz;);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    massfrac[CH4_ID] = prob_parm.Y_CH4_chamber;
    massfrac[O2_ID] = prob_parm.Y_O2_chamber;
    massfrac[N2_ID] = 1.0 - massfrac[O2_ID] - massfrac[CH4_ID];

    state(i,j,k,TEMP) = prob_parm.T_mean;

    // Read in turbulent fluctuations from MarcHDF
    state(i,j,k,VELX) = 0.0;
    state(i,j,k,VELY) = 0.0;
    state(i,j,k,VELZ) = 0.0;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i,j,k,TEMP), rho_cgs);
    state(i,j,k,DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i,j,k,TEMP), massfrac, h_cgs);
    state(i,j,k,RHOH) = h_cgs * 1.0e-4 * state(i,j,k,DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
       state(i,j,k,FIRSTSPEC+n) = massfrac[n] * state(i,j,k,DENSITY);
    }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
time_eta(
  amrex::Real time, amrex::Real start, amrex::Real duration, amrex::Real tau)
{
  amrex::Real eta_time;
  eta_time = 0.5 * (std::tanh((time - start) / (duration * tau)) *
                    std::tanh((-time + start + duration) / (duration * tau)) +
                    1.0);
  eta_time = std::min(std::max(eta_time, 0.0), 1.0);
  return eta_time;
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
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
  auto eos = pele::physics::PhysicsType::eos();

  // Chamber composition
  amrex::Real massfrac_ch[NUM_SPECIES] = {0.0};
  massfrac_ch[CH4_ID] = prob_parm.Y_CH4_chamber;
  massfrac_ch[O2_ID] = prob_parm.Y_O2_chamber;
  massfrac_ch[N2_ID] = 1.0 - massfrac_ch[O2_ID] - massfrac_ch[CH4_ID];

  // Pure fuel composition
  amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
  massfrac_fuel[NC12H26_ID] = 1.0;

  // Jet composition
  amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac_jet[n] = prob_parm.Z * massfrac_fuel[n] + (1.0 - prob_parm.Z) * massfrac_ch[n];
  }

  // Fill defaults: no slip wall, sort of.
  amrex::Real p_cgs = prob_parm.P_mean * 10.0;
  amrex::Real rho_cgs = 0.0;
  amrex::Real RhoH_temp = 0.0;
  s_ext[TEMP] = prob_parm.T_mean;
  eos.PYT2R(p_cgs, massfrac_ch, s_ext[TEMP], rho_cgs);
  s_ext[DENSITY] = rho_cgs * 1.0e3;
  eos.TY2H(s_ext[TEMP], massfrac_ch, RhoH_temp);
  s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC+n] = massfrac_ch[n] * s_ext[DENSITY];
  }

// Assuming this boundary is the top XY plane
// and in-plane angle is measured from X axis
  int isInjCell = 0;
  const amrex::Real dphi_angle = 2.0 * Pi / prob_parm.nholes;

  for (int nh = 1; nh <= prob_parm.nholes; nh++) {
    // Current hole center
    const amrex::Real phi_angle = nh * dphi_angle;
    const amrex::Real cx = prob_parm.centx + prob_parm.r_circ * cos(phi_angle);
    const amrex::Real cy = prob_parm.centy + prob_parm.r_circ * sin(phi_angle);

    // Injecting on an ellipse with R_minor = prob_parm.r_hole and R_major = prob_parm.r_hole / cos(prob_parm.cone_angle)
    // Local coordinates: 
    amrex::Real n_loc = (x[0] - cx) * cos(phi_angle) - (x[1] - cy) * sin(phi_angle);
    amrex::Real t_loc = (x[0] - cx) * sin(phi_angle) + (x[1] - cy) * cos(phi_angle);

    // Ellipse: major in normal, minor is tangential
    amrex::Real elrSq =  n_loc * n_loc / (  prob_parm.r_hole/cos(prob_parm.cone_angle* Pi / 180.0)
                                          * prob_parm.r_hole/cos(prob_parm.cone_angle* Pi / 180.0) )
                       + t_loc * t_loc / ( prob_parm.r_hole * prob_parm.r_hole );

    // We are using a tanH profile, past 1.2 keep no-slip walls conditions from above
    if ( elrSq < 1.2 ) {
      // Switch ON flag for injection face
      isInjCell = 1;

      // Tanh profile, assuming 10% tanh width
      amrex::Real eta = 0.5 * ( 1.0 - std::tanh((elrSq-1.0)/0.1));

      eta *= time_eta(time, prob_parm.inj_start, prob_parm.inj_dur, prob_parm.tau);

      // Composition
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++)
      {
        massfrac[n] = eta * massfrac_jet[n]+ (1.0-eta) * massfrac_ch[n];
      }

      // Velocity
      amrex::Real vx_in = 0.0, vy_in = 0.0, vz_in = 0.0;
      const amrex::Real vjet_pipeflow = prob_parm.vel_jet * eta;

      // find velocity vector
      vz_in =
        -vjet_pipeflow *
        cos(prob_parm.cone_angle * Pi / 180.0); // top XY plane
      const amrex::Real vr_in =
        vjet_pipeflow * sin(prob_parm.cone_angle * Pi / 180.0);
      vx_in = vr_in * cos(phi_angle);
      vy_in = vr_in * sin(phi_angle);

      s_ext[TEMP] = prob_parm.T_jet;
      eos.PYT2R(p_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;
      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
      // Add jet velocity to turbulence data
      s_ext[VELX] *= eta;
      s_ext[VELY] *= eta;
      s_ext[VELZ] *= eta;
      s_ext[VELX] += vx_in;
      s_ext[VELY] += vy_in;
      s_ext[VELZ] += vz_in;
    }
  }

  // If not a cell located on an injection, set velocity fluctuations to zero
  if (!isInjCell) {
     s_ext[VELX] = 0.0;
     s_ext[VELY] = 0.0;
     s_ext[VELZ] = 0.0;                                                                                                                                                      
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
zero_visc (int i, int j, int k,
           amrex::Array4<amrex::Real> const& beta,
           amrex::GeometryData const& geomdata,
           amrex::Box const& domainBox,
           const int  dir,
           const int  beta_comp,
           const int  nComp)
{
     amrex::ignore_unused(i,j,k,beta,geomdata,domainBox,dir,beta_comp,nComp);
    // We treat species when beta_comp == 0 and nComp == NUM_SPECIES
    // otherwise this routine could be called for other face diffusivity (Temp, velocity, ...)
}
#endif
