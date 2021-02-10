#ifndef PeleLM_K_H_
#define PeleLM_K_H_

#include <mechanism.h>
#include <EOS.H>
#include <Transport.H>
#include <TransportParams.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getTransportCoeff(int i, int j, int k,
                  amrex::Array4<const amrex::Real> const& rhoY,
                  amrex::Array4<const amrex::Real> const& T,
                  amrex::Array4<      amrex::Real> const& rhoDi,
                  amrex::Array4<      amrex::Real> const& lambda,
                  amrex::Array4<      amrex::Real> const& mu, 
                  TransParm const* trans_parm) noexcept
{
   using namespace amrex::literals;

   amrex::Real mwtinv[NUM_SPECIES] = {0.0};
   EOS::inv_molecular_weight(mwtinv);

   // Get rho & Y from rhoY
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }   
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }   

   amrex::Real Wbar = 0.0_rt;
   EOS::Y2WBAR(y, Wbar);

   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real rhoDi_cgs[NUM_SPECIES] = {0.0};
   amrex::Real lambda_cgs = 0.0_rt;
   amrex::Real mu_cgs = 0.0_rt;
   amrex::Real dummy_xi = 0.0_rt;
   amrex::Real Tloc = T(i,j,k);

   bool get_xi = false;
   bool get_mu = true;
   bool get_lam = true;
   bool get_Ddiag = true;
   transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
             rho, y, rhoDi_cgs, mu_cgs, dummy_xi, lambda_cgs, trans_parm);

   // Do CGS -> MKS conversions
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoDi(i,j,k,n) = rhoDi_cgs[n] * Wbar * mwtinv[n] * 1.0e-1_rt;
   }
   lambda(i,j,k) = lambda_cgs * 1.0e-5_rt;
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

#endif
