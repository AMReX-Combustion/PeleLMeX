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
                     amrex::Array4<amrex::Real> const& vel,
                     amrex::Array4<amrex::Real> const& rho,
                     amrex::Array4<amrex::Real> const& rhoY,
                     amrex::Array4<amrex::Real> const& rhoH,
                     amrex::Array4<amrex::Real> const& temp,
                     amrex::Array4<amrex::Real> const& aux,
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
    massfrac[O2_ID] = 0.233;
    massfrac[N2_ID] = 0.767;

    temp(i,j,k) = prob_parm.T_mean;

    amrex::Real r_xy = std::sqrt(x*x+y*y);
    amrex::Real angle = std::atan2(y,x);
    amrex::Real angular_vel = r_xy * 400.0 * std::tanh(300*(0.012-r_xy));

    vel(i,j,k,0) = angular_vel * std::sin(angle);
    vel(i,j,k,1) = -angular_vel * std::cos(angle);
    vel(i,j,k,2) = 0.0;

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, temp(i,j,k), rho_cgs);
    rho(i,j,k) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(temp(i,j,k), massfrac, h_cgs);
    rhoH(i,j,k) = h_cgs * 1.0e-4 * rho(i,j,k);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
       rhoY(i,j,k,n) = massfrac[n] * rho(i,j,k);
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
  pele::physics::PMF::PmfData::DataContainer const *pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
  auto eos = pele::physics::PhysicsType::eos();

  // Fill defaults: no slip wall, sort of.
  amrex::Real p_cgs = prob_parm.P_mean * 10.0;
  amrex::Real rho_cgs = 0.0;
  amrex::Real RhoH_temp = 0.0;
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[O2_ID] = 0.233;
  massfrac[N2_ID] = 0.767;
  s_ext[TEMP] = prob_parm.T_mean;
  eos.PYT2R(p_cgs, massfrac, s_ext[TEMP], rho_cgs);
  s_ext[DENSITY] = rho_cgs * 1.0e3;
  eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
  s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion
  for (int n = 0; n < NUM_SPECIES; n++) {
    s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
  }
  s_ext[VELX] = 0.0;
  s_ext[VELY] = 0.0;
  s_ext[VELZ] = 0.0;

// Assuming this boundary is the top XY plane
// and in-plane angle is measured from X axis
  const amrex::Real dphi_angle = 2.0 * Pi / prob_parm.nholes;

  for (int nh = 1; nh <= prob_parm.nholes; nh++) {
    const amrex::Real phi_angle = nh * dphi_angle;
    const amrex::Real hole_cx = prob_parm.centx + prob_parm.r_circ * cos(phi_angle);
    const amrex::Real hole_cy = prob_parm.centy + prob_parm.r_circ * sin(phi_angle);
    if ( ((x[0] - hole_cx) * (x[0] - hole_cx) +
          (x[1] - hole_cy) * (x[1] - hole_cy)) <
         (prob_parm.r_hole * prob_parm.r_hole)) {

      const amrex::Real rSq = (x[0] - hole_cx) * (x[0] - hole_cx) +
                              (x[1] - hole_cy) * (x[1] - hole_cy);
      const amrex::Real rad = sqrt(rSq);


      // subsonic velocity inlet
      amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
      massfrac_jet[NC7H16_ID] = 1.0;

      amrex::Real vx_in = 0.0, vy_in = 0.0, vz_in = 0.0;
      const amrex::Real vjet_pipeflow =
          prob_parm.vel_jet *
          (1.0 - rSq / (prob_parm.r_hole * prob_parm.r_hole));

      // find velocity vector
      vz_in =
        -vjet_pipeflow *
        cos(prob_parm.cone_angle * Pi / 180.0); // top XY plane
      const amrex::Real vr_in =
        vjet_pipeflow * sin(prob_parm.cone_angle * Pi / 180.0);
      vx_in = vr_in * cos(phi_angle);
      vy_in = vr_in * sin(phi_angle);

      s_ext[TEMP] = prob_parm.T_jet;
      eos.PYT2R(p_cgs, massfrac_jet, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;
      eos.TY2H(s_ext[TEMP], massfrac_jet, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion
      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC+n] = massfrac_jet[n] * s_ext[DENSITY];
      }
      s_ext[VELX] = vx_in;
      s_ext[VELY] = vy_in;
      s_ext[VELZ] = vz_in;
    }
  }
}
#endif
