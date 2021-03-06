#ifndef PELELM_PROB_H_
#define PELELM_PROB_H_

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLM_Index.H>
#include <pelelm_prob_parm.H>
#include <pmf.H>
#include <pmf_data.H>
#include <PelePhysics.H>
#include <EOS_Extension.H>
#include <PeleLMEF_Constants.H>

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
                     amrex::Array4<amrex::Real> const& nE,
                     amrex::Array4<amrex::Real> const& phiV,
                     amrex::GeometryData const& geomdata,
                     ProbParm const& prob_parm,
                     PmfData const * pmf_data)
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
    amrex::Real molefrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real pert;
    amrex::Real y1, y2;

    pert = 0.0;
    if (prob_parm.pertmag > 0.0) 
    {
       pert = prob_parm.pertmag * 
             (1.0 * std::sin(2 * Pi * 4 * x / Lx) +
              1.023 * std::sin(2 * Pi * 2 * (x - 0.004598) / Lx) +
              0.945 * std::sin(2 * Pi * 3 * (x - 0.00712435) / Lx) +
              1.017 * std::sin(2 * Pi * 5 * (x - 0.0033) / Lx) +
              0.982 * std::sin(2 * Pi * 5 * (x - 0.014234) / Lx));
    }

    y1 = (y - prob_parm.standoff - 0.5*dx[1] + pert);
    y2 = (y - prob_parm.standoff + 0.5*dx[1] + pert);

    y1 += 0.5*dx[1];
    PMF::pmf(pmf_data,y1, y1, pmf_vals);

    temp(i,j,k) = pmf_vals[0];;

    for (int n = 0; n < NUM_SPECIES; n++){
      massfrac[n] = pmf_vals[3 + n];
    }

    amrex::Real mw[NUM_SPECIES] = {0.0};
    eos.molecular_weight(mw);
    amrex::Real zk[NUM_SPECIES] = {0.0};
    pele::physics::eos::charge_mass(zk);
    for (int n = 0; n < NUM_SPECIES; n++) {
	    zk[n] *= 1000.0;    // CGS->MKS
	 }
    massfrac[E_ID] = 0.0;
    amrex::Real sum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = std::max(0.0,massfrac[n]);
        if ( n != N2_ID ) sum += massfrac[n];
    }
    massfrac[N2_ID] = 1.0 - sum;
    amrex::Real chargeDist = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
       chargeDist += massfrac[n] * zk[n];
    }
    if ( chargeDist < 0.0 ) {
       amrex::Real sum_pos_ion = 0.0;
       amrex::Real mean_zk_pos_ion = 0.0;
       for (int n = 0; n < NUM_SPECIES; n++) {
          if ( zk[n] > 0.0 ) {
             sum_pos_ion += massfrac[n];
             mean_zk_pos_ion += massfrac[n] * zk[n];
          }
       }
       for (int n = 0; n < NUM_SPECIES; n++) {
          if ( zk[n] > 0.0 ) {
             massfrac[N2_ID] += massfrac[n] / sum_pos_ion * chargeDist / mean_zk_pos_ion * sum_pos_ion;
             massfrac[n] -= massfrac[n] / sum_pos_ion * chargeDist / mean_zk_pos_ion * sum_pos_ion;   
          }
       }

    }
    chargeDist = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
       chargeDist += massfrac[n] * zk[n];
    }
    nE(i,j,k) = std::max(1.0e-24,chargeDist/elemCharge);


    AMREX_D_TERM(vel(i,j,k,0) = 0.0;,
                 vel(i,j,k,1) = pmf_vals[1];,
                 vel(i,j,k,2) = 0.0);

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

    nE(i,j,k) *= rho(i,j,k);
    phiV(i,j,k) = 0.0;
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
  PmfData const *pmf_data)
{
   const amrex::Real* prob_lo = geomdata.ProbLo();
   amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
   amrex::Real molefrac[NUM_SPECIES] = {0.0};
   amrex::Real massfrac[NUM_SPECIES] = {0.0};

   auto eos = pele::physics::PhysicsType::eos();
   if (sgn == 1) {
      PMF::pmf(pmf_data,prob_lo[idir], prob_lo[idir], pmf_vals);

      AMREX_D_TERM(s_ext[VELX] = 0.0;,
                   s_ext[VELY] = pmf_vals[1];,
                   s_ext[VELZ] = 0.0);

      s_ext[TEMP] = pmf_vals[0];
      for (int n = 0; n < NUM_SPECIES; n++){
         massfrac[n] = pmf_vals[3 + n];
      }

      massfrac[N2_ID] += massfrac[E_ID] + massfrac[H3Op_ID] + massfrac[HCOp_ID] + massfrac[C2H3Op_ID];
      massfrac[E_ID] = 0.0;
      massfrac[H3Op_ID] = 0.0;
      massfrac[HCOp_ID] = 0.0;
      massfrac[C2H3Op_ID] = 0.0;

      amrex::Real rho_cgs, P_cgs, RhoH_temp;
      P_cgs = prob_parm.P_mean * 10.0;

      eos.PYT2R(P_cgs, massfrac, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];   // CGS -> MKS conversion

      for (int n = 0; n < NUM_SPECIES; n++) {
         s_ext[FIRSTSPEC+n] = massfrac[n] * s_ext[DENSITY];
      }
   } else { // Some charged species have inflow at the outlet
      for (int n = 0; n < NUM_SPECIES; n++) {
         s_ext[FIRSTSPEC+n] = 0.0;
      }
   }

   if (idir == 1) {  // y-dir
      if ( sgn == 1 ) { // lo
         s_ext[NE] = 0.0;
         s_ext[PHIV] = prob_parm.phiV_loy;
      } else {          // hi
         s_ext[NE] = 0.0;
         s_ext[PHIV] = prob_parm.phiV_hiy;
      }   
   }
}
#endif
