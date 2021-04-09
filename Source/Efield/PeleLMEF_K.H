#ifndef PeleLMEF_K_H_
#define PeleLMEF_K_H_

#include <PeleLMEF_Constants.H>
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
           amrex::Array4<      amrex::Real> const& rhoY,
           amrex::Array4<      amrex::Real> const& rhoD,
           amrex::Array4<const amrex::Real> const& T,
           amrex::Array4<      amrex::Real> const& Ksp) noexcept
{
   using namespace amrex::literals;
   amrex::Real rho = 0.0;
   for (int n = 0; n < NUM_SPECIES; n++) {
      rho += rhoY(i,j,k,n);
   }
   for (int n = 0; n < NUM_SPECIES; n++) {
      Ksp(i,j,k,n) = rhoD(i,j,k,n) * mwt[n] * 0.001 * a_zk[n] / ( rho * PP_RU_MKS * T(i,j,k) );
      // TODO if ( n == E_ID ) Ksp(i,j,k,n) = 0.0;         // Set the E species mobility to zero. Electron are treated separately.
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

#endif
