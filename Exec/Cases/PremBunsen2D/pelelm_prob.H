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

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
    amrex::Real u[3] = {0.0};
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real y1, y2;
    amrex::Real splitx = (prob_parm.is_sym) ? prob_lo[0]
                                            : prob_lo[0] + 0.5 * Lx;

    massfrac[N2_ID] = 0.767;
    massfrac[O2_ID] = 0.233;

    state(i,j,k,TEMP) = 298.0;

    if ( x >= splitx ) {
       if ( (x-splitx) < ( std::sqrt((y+0.0009) / 0.01))/100.0) {
          amrex::Real pl1_angle = std::atan(0.01/prob_parm.slot_width);
          y1 = - ( x - splitx - prob_parm.slot_width ) * std::cos(pl1_angle) + y * std::sin(pl1_angle);
          y2 =   ( x - splitx - prob_parm.slot_width ) * std::sin(pl1_angle) + y * std::cos(pl1_angle);
          y1 += prob_parm.standoff;
          y2 += prob_parm.standoff;

          pele::physics::PMF::pmf(pmf_data,y2, y2, pmf_vals);
          state(i,j,k,TEMP) = pmf_vals[0];
          for (int n = 0; n < NUM_SPECIES; n++){
            massfrac[n] = pmf_vals[3 + n];
          }
       }
    } else {
       if ( std::abs((x-splitx)) < ( std::sqrt((y+0.0009) / 0.01))/100.0) {
          amrex::Real pl1_angle = std::atan(-0.01/prob_parm.slot_width);
          y1 = ( x - splitx + prob_parm.slot_width ) * std::cos(pl1_angle) - y * std::sin(pl1_angle);
          y2 = ( x - splitx + prob_parm.slot_width ) * std::sin(pl1_angle) + y * std::cos(pl1_angle);
          y1 += prob_parm.standoff;
          y2 += prob_parm.standoff;

          pele::physics::PMF::pmf(pmf_data, y2, y2, pmf_vals);
          state(i,j,k,TEMP) = pmf_vals[0];
          for (int n = 0; n < NUM_SPECIES; n++){
            massfrac[n] = pmf_vals[3 + n];
          }
       }
    }

    AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                 state(i,j,k,VELY) = 0.1;,
                 state(i,j,k,VELZ) = 0.0);

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
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();
   const amrex::Real* prob_hi = geomdata.ProbHi();

   AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
                const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
                const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

   amrex::Real splitx = (prob_parm.is_sym) ? prob_lo[0]
                                           : prob_lo[0] + 0.5 * Lx;
   amrex::Real massfrac_air[NUM_SPECIES] = {0.0};
   amrex::Real massfrac_jet[NUM_SPECIES] = {0.0};
   massfrac_air[N2_ID] = 0.767;
   massfrac_air[O2_ID] = 0.233;
   massfrac_jet[N2_ID] = 0.725;
   massfrac_jet[O2_ID] = 0.22;
   massfrac_jet[CH4_ID] = 0.055;

   auto eos = pele::physics::PhysicsType::eos();

   for (int n = 0; n < NVAR; n++){
      s_ext[n] = 0.0;
   }

// Get zone
   int zone = -1;
   if ( std::abs(x[0]-splitx) < prob_parm.slot_width ) {
      zone = 3; // jet
   } else {
      zone = 1; // coflow
   }

   if (idir == 1 && sgn == 1) {
      s_ext[VELX] = 0.0;
      amrex::Real rad = 0.0;
      if ( zone == 3 ) {
         rad = std::abs(x[0]-splitx);
         s_ext[VELY] = prob_parm.Vin * ( 1.0 - (rad/prob_parm.slot_width) * (rad/prob_parm.slot_width) );
      } else if ( zone == 1 ) {
         rad = std::abs(x[0]-splitx);
         s_ext[VELY] = prob_parm.Vcoflow * ( 1.0 + (rad)/0.025);
         if ( rad <= prob_parm.slot_width*1.4 ) s_ext[VELY] = 0.0;
      }

      s_ext[TEMP] = 298.0;

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;
      if ( zone == 1 ) {
         eos.PYT2R(P_cgs, massfrac_air, s_ext[TEMP], rho_cgs);
         s_ext[DENSITY] = rho_cgs * 1.0e3;

         eos.TY2H(s_ext[TEMP], massfrac_air, RhoH_temp);
         s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

         for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[FIRSTSPEC+n] = massfrac_air[n] * s_ext[DENSITY];
         }
      } else if ( zone == 3 ) {
         eos.PYT2R(P_cgs, massfrac_jet, s_ext[TEMP], rho_cgs);
         s_ext[DENSITY] = rho_cgs * 1.0e3;

         eos.TY2H(s_ext[TEMP], massfrac_jet, RhoH_temp);
         s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

         for (int n = 0; n < NUM_SPECIES; n++) {
            s_ext[FIRSTSPEC+n] = massfrac_jet[n] * s_ext[DENSITY];
         }
      }
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
