#ifndef PeleLM_K_H_
#define PeleLM_K_H_

#include <mechanism.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getTransportCoeff(int i, int j, int k,
                  amrex::Array4<const amrex::Real> const& rhoY,
                  amrex::Array4<const amrex::Real> const& T,
                  amrex::Array4<      amrex::Real> const& rhoDi,
                  amrex::Array4<      amrex::Real> const& lambda,
                  amrex::Array4<      amrex::Real> const& mu,
                  pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                      pele::physics::PhysicsType::transport_type> const* trans_parm) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real mwtinv[NUM_SPECIES] = {0.0};
   eos.inv_molecular_weight(mwtinv);

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
   eos.Y2WBAR(y, Wbar);

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
   auto trans = pele::physics::PhysicsType::transport();
   trans.transport(get_xi, get_mu, get_lam, get_Ddiag, Tloc,
                   rho, y, rhoDi_cgs, mu_cgs, dummy_xi, lambda_cgs, trans_parm);

   // Do CGS -> MKS conversions
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoDi(i,j,k,n) = rhoDi_cgs[n] * Wbar * mwtinv[n] * 1.0e-1_rt;
   }
   lambda(i,j,k) = lambda_cgs * 1.0e-5_rt;
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getVelViscosity(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& mu,
                pele::physics::transport::TransParm<pele::physics::PhysicsType::eos_type,
                                                    pele::physics::PhysicsType::transport_type> const* trans_parm) noexcept
{
   using namespace amrex::literals;

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

   rho *= 1.0e-3_rt;                          // MKS -> CGS conversion
   amrex::Real temp = T(i,j,k);
   amrex::Real dummy_rhoDi[NUM_SPECIES] = {0.0};
   amrex::Real dummy_lambda = 0.0_rt;
   amrex::Real mu_cgs = 0.0_rt;
   amrex::Real dummy_xi = 0.0_rt;

   bool get_xi = false;
   bool get_mu = true;
   bool get_lam = false;
   bool get_Ddiag = false;
   auto trans = pele::physics::PhysicsType::transport();
   trans.transport(get_xi, get_mu, get_lam, get_Ddiag, temp,
                   rho, y, dummy_rhoDi, mu_cgs, dummy_xi, dummy_lambda, trans_parm);

   // Do CGS -> MKS conversions
   mu(i,j,k) = mu_cgs * 1.0e-1_rt;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getPGivenRTY(int i, int j, int k,
             amrex::Array4<const amrex::Real> const& rho,
             amrex::Array4<const amrex::Real> const& rhoY,
             amrex::Array4<const amrex::Real> const& T,
             amrex::Array4<      amrex::Real> const& P) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real rho_cgs = rho(i,j,k) * 0.001_rt;

   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   eos.RTY2P(rho_cgs, T(i,j,k), y, P(i,j,k));
   P(i,j,k) = P(i,j,k) * 0.1_rt;                         // CGS -> MKS conversion
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
compute_divu(int i, int j, int k,
             amrex::Array4<const amrex::Real> const& rhoY,
             amrex::Array4<const amrex::Real> const& T,
             amrex::Array4<const amrex::Real> const& specDiff,
             amrex::Array4<const amrex::Real> const& tempDiff,
             amrex::Array4<const amrex::Real> const& specEnthDiff,
             amrex::Array4<const amrex::Real> const& rhoYdot,
             amrex::Array4<const amrex::Real> const& extRhoY,
             amrex::Array4<const amrex::Real> const& extRhoH,
             amrex::Array4<      amrex::Real> const& divu,
             int do_react) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real mwtinv[NUM_SPECIES] = {0.0};
   eos.inv_molecular_weight(mwtinv);

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

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   amrex::Real Wbar = 0.0_rt;
   eos.Y2WBAR(y, Wbar);
   amrex::Real hi[NUM_SPECIES] = {0.0};
   eos.T2Hi(T(i,j,k), hi);
   cpmix *= 1.0e-4_rt;                                // CGS -> MKS conversion
   for (int n = 0; n < NUM_SPECIES; n++) {
      hi[n] *= 1.0e-4_rt;                             // CGS -> MKS conversion
   }

   amrex::Real denominv = 1.0_rt / ( rho * cpmix * T(i,j,k) );
   divu(i,j,k) = ( specEnthDiff(i,j,k) + tempDiff(i,j,k) + extRhoH(i,j,k) ) * denominv;
   for (int n = 0; n < NUM_SPECIES; n++) {
      amrex::Real specTerm = specDiff(i,j,k,n) + extRhoY(i,j,k,n);
      if (do_react) specTerm += rhoYdot(i,j,k,n);
      divu(i,j,k) +=  specTerm * ( mwtinv[n] * Wbar * rhoinv - hi[n] * denominv );
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
extFluxDivergence_K(int i, int j, int k,
                    int ncomp,
                    AMREX_D_DECL(amrex::Array4<const amrex::Real> const& fluxX,
                                 amrex::Array4<const amrex::Real> const& fluxY,
                                 amrex::Array4<const amrex::Real> const& fluxZ),
                    amrex::Array4<const amrex::Real>              const& vol,
                    const amrex::Real                                  & scaling,
                    amrex::Array4<      amrex::Real>              const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) = AMREX_D_TERM (  ( fluxX(i+1,j,k,n) - fluxX(i,j,k,n) ),
                                   + ( fluxY(i,j+1,k,n) - fluxY(i,j,k,n) ),
                                   + ( fluxZ(i,j,k+1,n) - fluxZ(i,j,k,n) ));
      div(i,j,k,n) *= factor;
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
EB_intFluxDivergence_K(int i, int j, int k,
                       int ncomp,
                       AMREX_D_DECL(amrex::Array4<const amrex::Real> const& fluxX,
                                    amrex::Array4<const amrex::Real> const& fluxY,
                                    amrex::Array4<const amrex::Real> const& fluxZ),
                       AMREX_D_DECL(amrex::Array4<const amrex::Real> const& afracX,
                                    amrex::Array4<const amrex::Real> const& afracY,
                                    amrex::Array4<const amrex::Real> const& afracZ),
                       AMREX_D_DECL(amrex::Real areax,
                                    amrex::Real areay,
                                    amrex::Real areaz),
                       amrex::Array4<const amrex::Real>              const& vol,
                       const amrex::Real                                  & scaling,
                       amrex::Array4<      amrex::Real>              const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) = AMREX_D_TERM (  (   areax * afracX(i+1,j,k) * fluxX(i+1,j,k,n)
                                       - areax * afracX(i  ,j,k) * fluxX(i  ,j,k,n) ),
                                   + (   areay * afracY(i,j+1,k) * fluxY(i,j+1,k,n)
                                       - areay * afracY(i,j  ,k) * fluxY(i,j  ,k,n) ),
                                   + (   areaz * afracZ(i,j,k+1) * fluxZ(i,j,k+1,n)
                                       - areaz * afracZ(i,j,k  ) * fluxZ(i,j,k  ,n) ));
      div(i,j,k,n) *= factor;
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
intFluxDivergence_K(int i, int j, int k,
                    int ncomp,
                    AMREX_D_DECL(amrex::Array4<const amrex::Real> const& fluxX,
                                 amrex::Array4<const amrex::Real> const& fluxY,
                                 amrex::Array4<const amrex::Real> const& fluxZ),
                    AMREX_D_DECL(amrex::Real areax,
                                 amrex::Real areay,
                                 amrex::Real areaz),
                    amrex::Array4<const amrex::Real>              const& vol,
                    const amrex::Real                                  & scaling,
                    amrex::Array4<      amrex::Real>              const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) = AMREX_D_TERM (  (   areax * fluxX(i+1,j,k,n)
                                       - areax * fluxX(i  ,j,k,n) ),
                                   + (   areay * fluxY(i,j+1,k,n)
                                       - areay * fluxY(i,j  ,k,n) ),
                                   + (   areaz * fluxZ(i,j,k+1,n)
                                       - areaz * fluxZ(i,j,k  ,n) ));
      div(i,j,k,n) *= factor;
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
intFluxDivergence_rz_K(int i, int j, int k,
                       int ncomp,
                       AMREX_D_DECL(amrex::Array4<const amrex::Real> const& fluxX,
                                    amrex::Array4<const amrex::Real> const& fluxY,
                                    amrex::Array4<const amrex::Real> const& fluxZ),
                       amrex::Array4<const amrex::Real>              const& ax,
                       amrex::Array4<const amrex::Real>              const& ay,
                       amrex::Array4<const amrex::Real>              const& vol,
                       const amrex::Real                                  & scaling,
                       amrex::Array4<      amrex::Real>              const& div) noexcept
{
   using namespace amrex::literals;

   amrex::Real factor = scaling / vol(i,j,k);

   for (int n = 0; n < ncomp; n++) {
      div(i,j,k,n) =  (  ax(i+1,j,k) * fluxX(i+1,j,k,n)
                       - ax(i  ,j,k) * fluxX(i  ,j,k,n))
                    + (  ay(i,j+1,k) * fluxY(i,j+1,k,n)
                       - ay(i,j  ,k) * fluxY(i,j  ,k,n));
      div(i,j,k,n) *= factor;
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getMwmixGivenRY(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rho,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<      amrex::Real> const& Mwmix) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   eos.Y2WBAR(y, Mwmix(i,j,k));
   Mwmix(i,j,k) = Mwmix(i,j,k) * 0.001_rt;            // CGS -> MKS conversion
}


AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
cen2edg_cpp(int i, int j, int k,
            int dir,
            int ncomp,
            const bool                       do_harmonic,
            const bool                       on_lo,
            const bool                       on_hi,
            amrex::Array4<const amrex::Real> const& cfab,
            amrex::Array4<      amrex::Real> const& efab) noexcept
{
   using namespace amrex::literals;

   // Default behavior: id_l -1 in dir
   // Handle the BCs
   // Need -1 in id_l and id_h in dir on low Dirichlet BC
   // Need nothing on high Dirichlet BC
   int  id_l[3] = {i,j,k};
   int  id_h[3] = {i,j,k};
   if (!on_lo ) {
      if ( !on_hi ) {
         id_l[dir] -= 1;
      }
   } else {
      id_l[dir] -= 1;
      id_h[dir] -= 1;
   }
   if (do_harmonic) {
      for (int n = 0; n < ncomp; n++) {
         if ( ( cfab(id_l[0],id_l[1],id_l[2],n) * cfab(id_h[0],id_h[1],id_h[2],n) ) > 0.0_rt ) {
            efab(i,j,k,n) = 2.0_rt * ( cfab(id_l[0],id_l[1],id_l[2],n) * cfab(id_h[0],id_h[1],id_h[2],n) )
                                   / ( cfab(id_l[0],id_l[1],id_l[2],n) + cfab(id_h[0],id_h[1],id_h[2],n) );
         } else {
            efab(i,j,k,n) = 0.0_rt;
         }
      }
   } else {
      for (int n = 0; n < ncomp; n++) {
         efab(i,j,k,n) = 0.5_rt*(cfab(id_l[0],id_l[1],id_l[2],n) + cfab(id_h[0],id_h[1],id_h[2],n));
      }
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
repair_flux(int i, int j, int k,
            int dir,
            const bool                       on_lo,
            const bool                       on_hi,
            amrex::Array4<const amrex::Real> const& rhoY,
            amrex::Array4<      amrex::Real> const& flux ) noexcept
{
   using namespace amrex::literals;

   // Handle the BCs : need a -1 on id_l in dir away from the BC
   // Need -1 in id_l and id_h in dir on low Dirichlet BC
   // Need nothing on high Dirichlet BC
   int  id_l[3] = {i,j,k};
   int  id_h[3] = {i,j,k};
   if ( !on_lo ) {
      if ( !on_hi ) {
        id_l[dir] -= 1;
      }
   } else {
      id_l[dir] -= 1;
      id_h[dir] -= 1;
   }

   amrex::Real sumFlux = 0.0_rt;
   amrex::Real sumRhoYe = 0.0_rt;
   amrex::Real RhoYe[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      sumFlux += flux(i,j,k,n);
      RhoYe[n] = 0.5_rt * (rhoY(id_l[0],id_l[1],id_l[2],n) + rhoY(id_h[0],id_h[1],id_h[2],n));
      sumRhoYe += RhoYe[n];
   }
   sumRhoYe = 1.0_rt / sumRhoYe;
   for (int n = 0; n < NUM_SPECIES; n++) {
      flux(i,j,k,n) -= sumFlux * RhoYe[n] * sumRhoYe;
   }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
repair_flux_eb(int i, int j, int k,
               int dir,
               const bool                       on_lo,
               const bool                       on_hi,
               amrex::Array4<const amrex::Real> const& rhoY,
               amrex::Array4<const amrex::Real> const& rhoYed,
               amrex::Array4<const amrex::Real> const& areafrac,
               amrex::Array4<      amrex::Real> const& flux ) noexcept
{
   using namespace amrex::literals;

   // Handle the BCs : need a -1 on id_l in dir on Dirichlet BC
   int  id_l[3] = {i,j,k};
   id_l[dir] -= 1;
   if ( areafrac(i,j,k) > 0.0_rt ) {
      amrex::Real sumFlux = 0.0_rt;
      amrex::Real sumRhoYe = 0.0_rt;
      amrex::Real RhoYe[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
         sumFlux += flux(i,j,k,n);
         RhoYe[n] = rhoYed(i,j,k,n);
         if ( on_lo ) RhoYe[n] = rhoY(id_l[0],id_l[1],id_l[2],n);
         if ( on_hi ) RhoYe[n] = rhoY(i,j,k,n);
         sumRhoYe += RhoYe[n];
      }
      sumRhoYe = 1.0_rt / sumRhoYe;
      for (int n = 0; n < NUM_SPECIES; n++) {
         flux(i,j,k,n) -= sumFlux * RhoYe[n] * sumRhoYe;
      }
   } else {
      for (int n = 0; n < NUM_SPECIES; n++) {
         flux(i,j,k,n) = 0.0_rt;
      }
   }

}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getHGivenT(int i, int j, int k,
           amrex::Array4<const amrex::Real> const& T,
           amrex::Array4<      amrex::Real> const& Hi) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real hi_spec[NUM_SPECIES] = {0.0_rt};
   eos.T2Hi(T(i,j,k), hi_spec);
   for (int n = 0; n < NUM_SPECIES; n++) {
      Hi(i,j,k,n) = hi_spec[n] * 0.0001_rt;            // CGS -> MKS conversion
   }
   
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
makeVelForce(int i, int j, int k,
             int                                     is_incomp,
             const amrex::Real                     & rho_incomp,
             int                                     pseudo_gravity,
             const amrex::Real                     & /*time*/,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const gravity,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const gp0,
             const amrex::Real                     & dV_control,
             amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const /*dx*/,
             amrex::Array4<const amrex::Real> const& vel,
             amrex::Array4<const amrex::Real> const& rho,
             amrex::Array4<const amrex::Real> const& /*rhoY*/,
             amrex::Array4<const amrex::Real> const& /*rhoh*/,
             amrex::Array4<const amrex::Real> const& /*temp*/,
             amrex::Array4<const amrex::Real> const& extmom,
             amrex::Array4<const amrex::Real> const& extrho,
             amrex::Array4<      amrex::Real> const& force) noexcept
{
   using namespace amrex::literals;

   // Switch between incompressible/low-Mach rhos
   amrex::Real rho_lcl = 0.0;
   if ( is_incomp ) {
      rho_lcl = rho_incomp;
   } else {
      rho_lcl = rho(i,j,k);
   }

   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      force(i,j,k,idim) = -gp0[idim] + gravity[idim] * rho_lcl;
   }

   // TODO direction of pseudo grav should be runtime params
   if (pseudo_gravity) {
      force(i,j,k,AMREX_SPACEDIM-1) += dV_control * rho_lcl;
   }
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
     force(i,j,k,idim) += extmom(i,j,k,idim) - vel(i,j,k,idim) * extrho(i,j,k);
   }
}

// This kernel returns directly rhoHmix
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getRHmixGivenTY(int i, int j, int k,
                amrex::Array4<const amrex::Real> const& rho,
                amrex::Array4<const amrex::Real> const& rhoY,
                amrex::Array4<const amrex::Real> const& T,
                amrex::Array4<      amrex::Real> const& Hmix) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   eos.TY2H(T(i,j,k), y, Hmix(i,j,k));
   Hmix(i,j,k) = Hmix(i,j,k) * 1.0e-4_rt * rho(i,j,k);   // CGS -> MKS conversion
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getTfromHY(int i, int j, int k,
           amrex::Array4<const amrex::Real> const& rho,
           amrex::Array4<const amrex::Real> const& rhoY,
           amrex::Array4<const amrex::Real> const& rhoH,
           amrex::Array4<      amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real rhoinv = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }
   amrex::Real h_cgs = rhoH(i,j,k) * rhoinv * 1.0e4_rt;
   amrex::Real T_loc = T(i,j,k);
   eos.HY2T(h_cgs, y, T_loc);
   T(i,j,k) = T_loc;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
getCpmixGivenRYT(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rho,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& T,
                 amrex::Array4<      amrex::Real> const& cpmix) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real rhoinv  = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   eos.TY2Cp(T(i,j,k), y, cpmix(i,j,k));
   cpmix(i,j,k) *= 0.0001_rt;                         // CGS -> MKS conversion
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
buildAdvectionForcing(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& rho,
                      amrex::Array4<const amrex::Real> const& rhoY,
                      amrex::Array4<const amrex::Real> const& T,
                      amrex::Array4<const amrex::Real> const& dn,
                      amrex::Array4<const amrex::Real> const& ddn,
                      amrex::Array4<const amrex::Real> const& r,
                      amrex::Array4<const amrex::Real> const& extRhoY,
                      amrex::Array4<const amrex::Real> const& extRhoH,
                      amrex::Real                      const& dp0dt,
                      int                              const& closed_chamber,
                      int                                     do_react,
                      amrex::Array4<      amrex::Real> const& forceY,
                      amrex::Array4<      amrex::Real> const& forceT) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
   // Get species enthalpy
   amrex::Real hi_spec[NUM_SPECIES] = {0.0};
   eos.T2Hi(T(i,j,k), hi_spec);

   // For species mass fractions
   amrex::Real rhoinv = 1.0_rt / rho(i,j,k);
   amrex::Real y[NUM_SPECIES] = {0.0};

   forceT(i,j,k) = dn(i,j,k,NUM_SPECIES) + ddn(i,j,k);
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;                                 // get y
      forceY(i,j,k,n) = dn(i,j,k,n);                                 // forceY = Dnk
      if (do_react) forceY(i,j,k,n) += r(i,j,k,n);                   //         + Rk 
      forceT(i,j,k) -= forceY(i,j,k,n) * hi_spec[n] * 0.0001_rt;     // forceT -= Sum{ hk. (Rk + Dnk) } with CGS -> MKS conversion of hk
      forceY(i,j,k,n) += extRhoY(i,j,k,n);
   }

   forceT(i,j,k) += extRhoH(i,j,k);
   if ( closed_chamber == 1 ) forceT(i,j,k) += dp0dt;                // Add dp0/dt

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   amrex::Real cpmixinv = 1.0_rt / cpmix * 1.0e4_rt;                 // CGS -> MKS conversion
   forceT(i,j,k) *= rhoinv * cpmixinv;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
buildDiffusionForcing(int i, int j, int k,
                      amrex::Array4<const amrex::Real> const& dn,
                      amrex::Array4<const amrex::Real> const& ddn,
                      amrex::Array4<const amrex::Real> const& dnp1k,
                      amrex::Array4<const amrex::Real> const& ddnp1k,
                      amrex::Array4<const amrex::Real> const& r,
                      amrex::Array4<const amrex::Real> const& a,
                      amrex::Real                      const& dp0dt,
                      int                              const& closed_chamber,
                      int                                     do_react,
                      amrex::Array4<      amrex::Real> const& forceY,
                      amrex::Array4<      amrex::Real> const& forceT) noexcept
{
   using namespace amrex::literals;

   for (int n = 0; n < NUM_SPECIES; n++) {
      forceY(i,j,k,n) = a(i,j,k,n) + 0.5_rt * ( dn(i,j,k,n) - dnp1k(i,j,k,n) );
      if (do_react) forceY(i,j,k,n) += r(i,j,k,n);
   }
   forceT(i,j,k) = a(i,j,k,NUM_SPECIES) + 0.5_rt * (  dn(i,j,k,NUM_SPECIES) + ddn(i,j,k)
                             - dnp1k(i,j,k,NUM_SPECIES) - ddnp1k(i,j,k) );

   if ( closed_chamber == 1 ) forceT(i,j,k) += dp0dt;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
reactionRateRhoY(int i, int j, int k,
                 amrex::Array4<const amrex::Real> const& rhoY,
                 amrex::Array4<const amrex::Real> const& rhoH,
                 amrex::Array4<const amrex::Real> const& T,
                 amrex::Array4<      amrex::Real> const& rhoYdot) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();

   // Get rho & Ys from rhoYs.
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0_rt};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   // Get T from Y/H.
   amrex::Real Tloc = T(i,j,k);
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;      // Include MKS -> CGS conversion
   eos.HY2T(H, y, Tloc);
   
   // Get wdot.
   rho = rho * 0.001_rt;                                 // rho MKS -> CGS
   amrex::Real wdot[NUM_SPECIES] = {0.0_rt};
   eos.RTY2WDOT(rho, Tloc, y, wdot);
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoYdot(i,j,k,n) = wdot[n] * 1000.0_rt;            // CGS -> MKS conversion
   }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fabMinMax(int i, int j, int k,
          int           ncomp,
          amrex::Real   fmin,
          amrex::Real   fmax,
          amrex::Array4<amrex::Real> const& fab) noexcept
{
   using namespace amrex::literals;

   for (int n = 0; n < ncomp; n++) {
      fab(i,j,k,n) = amrex::max( fmin, amrex::min( fmax, fab(i,j,k,n) ) );
   }
}

AMREX_FORCE_INLINE
void
linearChemForcing(const amrex::Box& bx,
                  amrex::Array4<amrex::Real> const& rY_in,
                  amrex::Array4<amrex::Real> const& rY_src_in,
                  amrex::Array4<amrex::Real> const& T_in,
                  amrex::Array4<amrex::Real> const& rhoH_in,
                  amrex::Array4<amrex::Real> const& rhoH_src_in,
                  amrex::Array4<amrex::Real> const& FC_in,
                  amrex::Array4<amrex::Real> const& avgIR,
                  amrex::Real &dt_react)
{
   ParallelFor(bx,
   [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
   {
      // Integration of external rhoH source term
      rhoH_in(i,j,k,0) += dt_react * rhoH_src_in(i,j,k,0);

      // Integration of external rhoY source term and avgDown IR
      for (int n = 0; n < NUM_SPECIES; n++){
          rY_in(i,j,k,n) += dt_react * ( rY_src_in(i,j,k,n) + avgIR(i,j,k,n) );
          //TODO need this ? rY_in(i,j,k,n) = amrex::max(0.0,rY_in(i,j,k,n));
      }

      /* T update with rho, H and Y */
      // Get updated rho
      amrex::Real rho_loc = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
          rho_loc += rY_in(i,j,k,n);
      }
      // Get updated Ys
      amrex::Real Y_loc[NUM_SPECIES] = {0.0};
      for (int n = 0; n < NUM_SPECIES; n++) {
          Y_loc[n] = rY_in(i,j,k,n) / rho_loc;
      }
      // Get updated H
      amrex::Real H_loc = rhoH_in(i,j,k,0) / rho_loc;

      // Get curr estimate of T
      amrex::Real T_loc    = T_in(i,j,k,0);
      auto eos = pele::physics::PhysicsType::eos();
      eos.HY2T(H_loc,Y_loc,T_loc);
      T_in(i,j,k,0) = T_loc;

      /* Dummy */
      FC_in(i,j,k,0) = 0.0;
   });
}

// DivuCeiling type 1 est_divu_dt function
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
est_divu_dt_1(int i, int j, int k,
              amrex::Real                                  const& dtfactor,
              amrex::Real                                  const& rhomin,
              amrex::Array4<const amrex::Real>             const& rho,
              amrex::Array4<const amrex::Real>             const& divu) noexcept
{
   using namespace amrex::literals;

   amrex::Real dtcell = 1.0e20_rt;
   bool rho_is_okay = (rho(i,j,k) > rhomin);

   if ( divu(i,j,k) > 0.0_rt ) {
      dtcell = (rho_is_okay) ? ( 1.0_rt - rhomin / rho(i,j,k) ) / divu(i,j,k) : 1.0_rt / divu(i,j,k);
      dtcell *= dtfactor;
   }
   return dtcell;
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
getGammaInv(int i, int j, int k,
         amrex::Array4<const amrex::Real> const& rhoY,
         amrex::Array4<const amrex::Real> const& T) noexcept
{
   using namespace amrex::literals;

   auto eos = pele::physics::PhysicsType::eos();
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

   amrex::Real cpmix = 0.0_rt;
   eos.TY2Cp(T(i,j,k), y, cpmix);
   amrex::Real cvmix = 0.0_rt;
   eos.TY2Cv(T(i,j,k), y, cvmix);
   
   amrex::Real gammainv = cvmix / cpmix;

   return gammainv;
}
#endif
