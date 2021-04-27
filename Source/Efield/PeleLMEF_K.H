#ifndef PeleLMEF_K_H_
#define PeleLMEF_K_H_

#include <PeleLMEF_Constants.H>
#include <PelePhysics.H>
#include <mechanism.h>

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getKappaE(int i, int j, int k,
          amrex::Array4<amrex::Real> const& Ke) noexcept
{
      using namespace amrex::literals;

         Ke(i,j,k) = 0.4_rt;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getKappaSp(int i, int j, int k,
           amrex::Real const mwt[NUM_SPECIES],
           amrex::GpuArray<amrex::Real, NUM_SPECIES> a_zk,
           amrex::Array4<const amrex::Real> const& rhoY,
           amrex::Array4<      amrex::Real> const& rhoD,
           amrex::Array4<const amrex::Real> const& T,
           amrex::Array4<      amrex::Real> const& Ksp) noexcept
{
   using namespace amrex::literals;
   amrex::Real rho = 0.0;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   int ionCnt = 0;
   for (int n = 0; n < NUM_SPECIES; n++) {
      if ( a_zk[n] != 0.0 ) {
         Ksp(i,j,k,ionCnt) = rhoD(i,j,k,n) * mwt[n] * 0.001 * a_zk[n] / ( rho * PP_RU_MKS * T(i,j,k) );
         if ( n == E_ID ) Ksp(i,j,k,ionCnt) = 0.0;         // Set the E species mobility to zero. Electron are treated separately.
         ionCnt += 1;
      }
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getDiffE(int i, int j, int k,
         const amrex::Real &factor,
         amrex::Array4<const amrex::Real> const& T,
         amrex::Array4<const amrex::Real> const& Ke,
         amrex::Array4<      amrex::Real> const& De) noexcept
{
   using namespace amrex::literals;

   De(i,j,k) = Ke(i,j,k) * T(i,j,k) * factor;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
ef_edge_state_extdir(int i,int j,int k,
                     int dir,
                     const bool                       on_lo,
                     const bool                       on_hi,
                     amrex::Array4<amrex::Real const> const& state,
                     amrex::Array4<amrex::Real const> const& effvel ) noexcept
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

   amrex::Real state_ed;
   if ( effvel(i,j,k) > smallVel) {
      state_ed = state(id_l[0],id_l[1],id_l[2]);
   } else if ( effvel(i,j,k) < -smallVel) {
      state_ed = state(id_h[0],id_h[1],id_h[2]);
   } else {
      state_ed = 0.5_rt * ( state(id_l[0],id_l[1],id_l[2]) +
                            state(id_h[0],id_h[1],id_h[2]) );
   }
   return state_ed;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
cen2edg_upwind(int i, int j, int k,
               int dir,
               int ncomp,
               const bool                       on_lo,
               const bool                       on_hi,
               amrex::Array4<const amrex::Real> const& ueff,
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

   for (int n = 0; n < ncomp; n++) {
      if ( ueff(i,j,k) > smallVel ) {
         efab(i,j,k,n) = cfab(id_l[0],id_l[1],id_l[2],n);
      } else if ( ueff(i,j,k) < -smallVel ) {
         efab(i,j,k,n) = cfab(id_h[0],id_h[1],id_h[2],n);
      } else {
         efab(i,j,k,n) = 0.5_rt*(cfab(id_l[0],id_l[1],id_l[2],n) + cfab(id_h[0],id_h[1],id_h[2],n));
      }
   }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
reactionRateRhoY_EF(int i, int j, int k,
                    amrex::Array4<const amrex::Real> const& rhoY,
                    amrex::Array4<const amrex::Real> const& rhoH,
                    amrex::Array4<const amrex::Real> const& T,
                    amrex::Array4<const amrex::Real> const& nE,
                    amrex::Array4<const amrex::Real> const& mask,
                    amrex::Array4<      amrex::Real> const& rhoYdot,
                    amrex::Array4<      amrex::Real> const& nEdot) noexcept
{
   using namespace amrex::literals;

   // mask < 0 is EB covered cells: no reactions.
   if ( mask(i,j,k) < 0.0 ) {
      for (int n = 0; n < NUM_SPECIES; n++) {
         rhoYdot(i,j,k,n) = 0.0;
      }
      nEdot(i,j,k) = 0.0;
      return;
   }

   // Get rho from rhoY. Is there a better way ?
   amrex::Real rho = 0.0_rt;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   amrex::Real rhoinv = 1.0_rt / rho;
   amrex::Real y[NUM_SPECIES] = {0.0};
   for (int n = 0; n < NUM_SPECIES; n++) {
      y[n] = rhoY(i,j,k,n) * rhoinv;
   }

   // Get Y[e] from nE
   auto eos = pele::physics::PhysicsType::eos();
   amrex::Real mw[NUM_SPECIES] = {0.0};
   eos.molecular_weight(mw);
   y[E_ID] = nE(i,j,k) / Na * mw[E_ID] * 0.001 * rhoinv;

   // Get T from Y/H. Was done like that before.
   amrex::Real Tloc = T(i,j,k);
   amrex::Real H = rhoH(i,j,k) * rhoinv * 1.0e4_rt;      // Include MKS -> CGS conversion
   eos.HY2T(H, y, Tloc);

   // Get wdot. Don't like the temporary wdot there ...
   rho = rho * 0.001_rt;                                 // rho MKS -> CGS
   amrex::Real wdot[NUM_SPECIES] = {0.0};
   eos.RTY2WDOT(rho, Tloc, y, wdot);
   for (int n = 0; n < NUM_SPECIES; n++) {
      rhoYdot(i,j,k,n) = wdot[n] * 1000.0_rt;            // CGS -> MKS conversion
   }
   nEdot(i,j,k) = wdot[E_ID] * Na * 1.0e6;
   rhoYdot(i,j,k,E_ID) = 0.0;
}
#endif
