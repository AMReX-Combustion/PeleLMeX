#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
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
                     pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
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

    AMREX_D_TERM(const amrex::Real x_c = prob_lo[0] + 0.5*Lx;,
                 const amrex::Real y_c = prob_lo[1] + 0.5*Ly;,
                 const amrex::Real z_c = prob_lo[2] + 0.5*Lz);

    auto eos = pele::physics::PhysicsType::eos();
    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    amrex::Real massfrac[NUM_SPECIES] = {0.0};

    // Species IDs
#ifdef O2_ID
    int sp1_ID = O2_ID;
#endif
#ifdef N2a_ID
    int sp1_ID = N2a_ID;
#endif
    int sp2_ID = N2_ID;

    massfrac[sp1_ID] = 0.233;
    massfrac[sp2_ID] = 0.767;

    if ( prob_parm.probType == 0 ) { // CoVo
       const amrex::Real deltax = x - prob_parm.xvort;
       const amrex::Real deltay = y - prob_parm.yvort;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = prob_parm.rvort * prob_parm.rvort;
       const amrex::Real u_vort = -prob_parm.forcevort*deltay/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real v_vort = prob_parm.forcevort*deltax/r_sq * exp(-d_sq/r_sq/2.);
       const amrex::Real w_vort = 0.;

       switch(prob_parm.meanFlowDir) {
         case 1  :
            AMREX_D_TERM(state(i,j,k,VELX) = prob_parm.meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -1  :
            AMREX_D_TERM(state(i,j,k,VELX) = -prob_parm.meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 2  :
            AMREX_D_TERM(state(i,j,k,VELX) = u_vort;,
                         state(i,j,k,VELY) = prob_parm.meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -2  :
            AMREX_D_TERM(state(i,j,k,VELX) = u_vort;,
                         state(i,j,k,VELY) = -prob_parm.meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 3  :
            AMREX_D_TERM(state(i,j,k,VELX) = prob_parm.meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = prob_parm.meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -3  :
            AMREX_D_TERM(state(i,j,k,VELX) = -prob_parm.meanFlowMag + u_vort;,
                         state(i,j,k,VELY) = -prob_parm.meanFlowMag + v_vort;,
                         state(i,j,k,VELZ) = 0.0);
            break;
       }

       if (is_incompressible) return;

       state(i,j,k,TEMP) = prob_parm.T_mean;
    } else if ( prob_parm.probType == 1 ) { // CoGau

       switch(prob_parm.meanFlowDir) {
         case 1  :
            AMREX_D_TERM(state(i,j,k,VELX) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -1  :
            AMREX_D_TERM(state(i,j,k,VELX) = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 2  :
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -2  :
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 3  :
            AMREX_D_TERM(state(i,j,k,VELX) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -3  :
            AMREX_D_TERM(state(i,j,k,VELX) = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
       }

       if (is_incompressible) return;

       state(i,j,k,TEMP) = prob_parm.T_mean;
       const amrex::Real deltax = x - prob_parm.xgauss;
       const amrex::Real deltay = y - prob_parm.ygauss;
       const amrex::Real d_sq = deltax*deltax + deltay*deltay;
       const amrex::Real r_sq = prob_parm.rgauss * prob_parm.rgauss;
       if ( prob_parm.gauss_type == 0 ) { // Spec
          massfrac[sp1_ID] += massfrac[sp1_ID] * prob_parm.ampgauss * std::exp(-d_sq/r_sq);
          massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
       } else if ( prob_parm.gauss_type == 1 ) { // Temp
          state(i,j,k,TEMP) = prob_parm.T_mean * ( 1.0 + prob_parm.ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
    } else if ( prob_parm.probType == 3 ) { // CoTanH

       switch(prob_parm.meanFlowDir) {
         case 1  :
            AMREX_D_TERM(state(i,j,k,VELX) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -1  :
            AMREX_D_TERM(state(i,j,k,VELX) = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = 0.0;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 2  :
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -2  :
            AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                         state(i,j,k,VELY) = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case 3  :
            AMREX_D_TERM(state(i,j,k,VELX) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELY) = prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ) = 0.0);
            break;
         case -3  :
            AMREX_D_TERM(state(i,j,k,VELX)  = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELY)  = -prob_parm.meanFlowMag;,
                         state(i,j,k,VELZ)  = 0.0);
            break;
       }

       if (is_incompressible) return;

       state(i,j,k,TEMP) = prob_parm.T_mean;
       const amrex::Real deltax = x - prob_parm.xgauss;
       const amrex::Real deltay = y - prob_parm.ygauss;
       const amrex::Real rad = std::sqrt(deltax*deltax + deltay*deltay);
       const amrex::Real eta = 0.5*(1.0 - std::tanh(2.0*(rad-prob_parm.rgauss)/(0.05*prob_parm.rgauss)));
       if ( prob_parm.gauss_type == 0 ) { // Spec
          massfrac[sp1_ID] += massfrac[sp1_ID] * prob_parm.ampgauss * eta;
          massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
       } else if ( prob_parm.gauss_type == 1 ) { // Temp
          state(i,j,k,TEMP) = prob_parm.T_mean * ( 1.0 + prob_parm.ampgauss * eta );
       } else {
          amrex::Abort("Unknown tanh_type: should be either Temp or Spec");
       }
    } else if ( prob_parm.probType == 2 ) { // DiffGau

       AMREX_D_TERM(state(i,j,k,VELX) = 0.0;,
                    state(i,j,k,VELY) = 0.0;,
                    state(i,j,k,VELZ) = 0.0);

       if (is_incompressible) return;

       state(i,j,k,TEMP) = prob_parm.T_mean;
       // Set up the Gaussian as the solution of the diffusion of a delta dirac initial distribution
       // after time gaussTime of a diffusion process with constant diffusion coeff. gaussDiff. Use ampgauss
       // to avoid undershoot/overshoot of the 2-species mixture.
       const amrex::Real deltax = x - prob_parm.xgauss;
       const amrex::Real deltay = y - prob_parm.ygauss;
       amrex::Real denom_inv = prob_parm.ampgauss / std::sqrt(4*Pi*prob_parm.gaussTime*prob_parm.gaussDiff);
       const amrex::Real d_sq = deltax*deltax;
       const amrex::Real r_sq = 4.0*prob_parm.gaussTime*prob_parm.gaussDiff;
       if ( prob_parm.gauss_type == 0 ) { // Spec
          massfrac[sp1_ID] = denom_inv * std::exp(-d_sq/r_sq);
          massfrac[sp2_ID] = 1.0 - massfrac[sp1_ID];
       } else if ( prob_parm.gauss_type == 1 ) { //Temp
          // TODO: find a better way, one with an analytical solution
          state(i,j,k,TEMP) = prob_parm.T_mean * ( 1.0 + prob_parm.ampgauss * std::exp(-d_sq/r_sq) );
       } else {
          amrex::Abort("Unknown gauss_type: should be either Temp or Spec");
       }
    }

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
  const amrex::Real* /*x[AMREX_SPACEDIM]*/,
  const int /*m_nAux*/,
  amrex::Real* /*s_ext[NVAR]*/,
  const int /*idir*/,
  const int /*sgn*/,
  const amrex::Real /*time*/,
  amrex::GeometryData const& /*geomdata*/,
  ProbParm const& /*prob_parm*/,
  pele::physics::PMF::PmfData::DataContainer const * /*pmf_data*/)
{
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
