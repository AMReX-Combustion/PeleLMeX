#include <PeleLMEFDeriveFunc.H>
#include <PeleLMEF_Constants.H>
#include <PeleLM_Index.H>
#include <PelePhysics.H>
#include <EOS_Extension.H>
#include <mechanism.H>

using namespace amrex;

void pelelm_derchargedist(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                          const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                          const amrex::Geometry& /*geomdata*/,
                          amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &/*bcrec*/, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& rhoY = statefab.const_array(FIRSTSPEC);
   auto const& nE   = statefab.const_array(NE);
   auto       der = derfab.array(dcomp);

   amrex::GpuArray<amrex::Real,NUM_SPECIES> zk;
   pele::physics::eos::charge_mass(zk.arr);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      der(i,j,k) = - nE(i,j,k) * elemCharge;
      for (int n = 0; n < NUM_SPECIES; ++n ) {
         der(i,j,k) += zk[n] * 1000.0 * rhoY(i,j,k,n); // CGS->MKS convertion of zk
      }
   });
}

void pelelm_derefx(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                   const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &bcrec, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& phiV = statefab.const_array(PHIV);
   auto       der = derfab.array(dcomp);

   const auto dxinv = geomdata.InvCellSizeArray();
   const auto domain = geomdata.Domain();
   amrex::Real factor = -0.5*dxinv[0];
   
   const auto bc_lo = bcrec[PHIV].lo(0);
   const auto bc_hi = bcrec[PHIV].hi(0);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && i <= domain.smallEnd(0) );
      bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && i >= domain.bigEnd(0) );
      der(i,j,k) = factor * ( phiV(i+1,j,k) - phiV(i-1,j,k) );
      if ( on_lo ) der(i,j,k) = factor * ( phiV(i+1,j,k) + phiV(i,j,k) - 2.0 * phiV(i-1,j,k) ) ;
      if ( on_hi ) der(i,j,k) = factor * ( 2.0 * phiV(i+1,j,k) - phiV(i,j,k) - phiV(i-1,j,k) ) ;
   });
}

void pelelm_derLorentzx(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                        const amrex::Geometry& geomdata,
                        amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &bcrec, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& phiV = statefab.const_array(PHIV);
   auto const& rhoY = statefab.const_array(FIRSTSPEC);
   auto const& nE   = statefab.const_array(NE);
   auto       der = derfab.array(dcomp);

   amrex::GpuArray<amrex::Real,NUM_SPECIES> zk;
   pele::physics::eos::charge_mass(zk.arr);

   const auto dxinv = geomdata.InvCellSizeArray();
   const auto domain = geomdata.Domain();
   amrex::Real factor = -0.5*dxinv[0];
   
   const auto bc_lo = bcrec[PHIV].lo(0);
   const auto bc_hi = bcrec[PHIV].hi(0);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      // Get gradient of PhiV
      bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && i <= domain.smallEnd(0) );
      bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && i >= domain.bigEnd(0) );
      amrex::Real EFx = factor * ( phiV(i+1,j,k) - phiV(i-1,j,k) );
      if ( on_lo ) EFx = factor * ( phiV(i+1,j,k) + phiV(i,j,k) - 2.0 * phiV(i-1,j,k) ) ;
      if ( on_hi ) EFx = factor * ( 2.0 * phiV(i+1,j,k) - phiV(i,j,k) - phiV(i-1,j,k) ) ;

      // Assemble Lorentz force in X
      der(i,j,k) = - nE(i,j,k) * elemCharge * EFx;
      for (int n = 0; n < NUM_SPECIES; ++n ) {
         der(i,j,k) += zk[n] * 1000.0 * rhoY(i,j,k,n) * EFx; // CGS->MKS convertion of zk
      }
   });
}

#if (AMREX_SPACEDIM > 1)
void pelelm_derefy(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                   const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &bcrec, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& phiV = statefab.const_array(PHIV);
   auto       der = derfab.array(dcomp);

   const auto dxinv = geomdata.InvCellSizeArray();
   const auto domain = geomdata.Domain();
   amrex::Real factor = -0.5*dxinv[1];
   
   const auto bc_lo = bcrec[PHIV].lo(1);
   const auto bc_hi = bcrec[PHIV].hi(1);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && j <= domain.smallEnd(1) );
      bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && j >= domain.bigEnd(1) );
      der(i,j,k) = factor * ( phiV(i,j+1,k) - phiV(i,j-1,k) );
      if ( on_lo ) der(i,j,k) = factor * ( phiV(i,j+1,k) + phiV(i,j,k) - 2.0 * phiV(i,j-1,k) ) ;
      if ( on_hi ) der(i,j,k) = factor * ( 2.0 * phiV(i,j+1,k) - phiV(i,j,k) - phiV(i,j-1,k) ) ;
   });
}

void pelelm_derLorentzy(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                        const amrex::Geometry& geomdata,
                        amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &bcrec, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& phiV = statefab.const_array(PHIV);
   auto const& rhoY = statefab.const_array(FIRSTSPEC);
   auto const& nE   = statefab.const_array(NE);
   auto       der = derfab.array(dcomp);

   amrex::GpuArray<amrex::Real,NUM_SPECIES> zk;
   pele::physics::eos::charge_mass(zk.arr);

   const auto dxinv = geomdata.InvCellSizeArray();
   const auto domain = geomdata.Domain();
   amrex::Real factor = -0.5*dxinv[1];
   
   const auto bc_lo = bcrec[PHIV].lo(1);
   const auto bc_hi = bcrec[PHIV].hi(1);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      // Get gradient of PhiV
      bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && j <= domain.smallEnd(1) );
      bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && j >= domain.bigEnd(1) );
      amrex::Real EFy = factor * ( phiV(i,j+1,k) - phiV(i,j-1,k) );
      if ( on_lo ) EFy = factor * ( phiV(i,j+1,k) + phiV(i,j,k) - 2.0 * phiV(i,j-1,k) ) ;
      if ( on_hi ) EFy = factor * ( 2.0 * phiV(i,j-1,k) - phiV(i,j,k) - phiV(i,j-1,k) ) ;

      // Assemble Lorentz force in Y
      der(i,j,k) = - nE(i,j,k) * elemCharge * EFy;
      for (int n = 0; n < NUM_SPECIES; ++n ) {
         der(i,j,k) += zk[n] * 1000.0 * rhoY(i,j,k,n) * EFy; // CGS->MKS convertion of zk
      }
   });
}

#if (AMREX_SPACEDIM > 2)
void pelelm_derefz(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                   const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                   const amrex::Geometry& geomdata,
                   amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &bcrec, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& phiV = statefab.const_array(PHIV);
   auto       der = derfab.array(dcomp);

   const auto dxinv = geomdata.InvCellSizeArray();
   const auto domain = geomdata.Domain();
   amrex::Real factor = -0.5*dxinv[2];
   
   const auto bc_lo = bcrec[PHIV].lo(2);
   const auto bc_hi = bcrec[PHIV].hi(2);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && k <= domain.smallEnd(2) );
      bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && k >= domain.bigEnd(2) );
      der(i,j,k) = factor * ( phiV(i,j,k+1) - phiV(i,j,k-1) );
      if ( on_lo ) der(i,j,k) = factor * ( phiV(i,j,k+1) + phiV(i,j,k) - 2.0 * phiV(i,j,k-1) ) ;
      if ( on_hi ) der(i,j,k) = factor * ( 2.0 * phiV(i,j,k+1) - phiV(i,j,k) - phiV(i,j,k-1) ) ;
   });
}

void pelelm_derLorentzz(const amrex::Box& bx, amrex::FArrayBox& derfab, int dcomp, int /*ncomp*/,
                        const amrex::FArrayBox& statefab, const amrex::FArrayBox& /*pressfab*/,
                        const amrex::Geometry& geomdata,
                        amrex::Real /*time*/, const amrex::Vector<amrex::BCRec> &bcrec, int /*level*/)
{
   AMREX_ASSERT(derfab.box().contains(bx));
   AMREX_ASSERT(statefab.box().contains(bx));
   auto const& phiV = statefab.const_array(PHIV);
   auto const& rhoY = statefab.const_array(FIRSTSPEC);
   auto const& nE   = statefab.const_array(NE);
   auto       der = derfab.array(dcomp);

   amrex::GpuArray<amrex::Real,NUM_SPECIES> zk;
   pele::physics::eos::charge_mass(zk.arr);

   const auto dxinv = geomdata.InvCellSizeArray();
   const auto domain = geomdata.Domain();
   amrex::Real factor = -0.5*dxinv[2];
   
   const auto bc_lo = bcrec[PHIV].lo(2);
   const auto bc_hi = bcrec[PHIV].hi(2);
   
   amrex::ParallelFor(bx,
   [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
   {
      // Get gradient of PhiV
      bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && k <= domain.smallEnd(2) );
      bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && k >= domain.bigEnd(2) );
      amrex::Real EFz = factor * ( phiV(i,j,k+1) - phiV(i,j,k-1) );
      if ( on_lo ) EFz = factor * ( phiV(i,j,k+1) + phiV(i,j,k) - 2.0 * phiV(i,j,k-1) ) ;
      if ( on_hi ) EFz = factor * ( 2.0 * phiV(i,j,k+1) - phiV(i,j,k) - phiV(i,j,k-1) ) ;

      // Assemble Lorentz force in Z
      der(i,j,k) = - nE(i,j,k) * elemCharge * EFz;
      for (int n = 0; n < NUM_SPECIES; ++n ) {
         der(i,j,k) += zk[n] * 1000.0 * rhoY(i,j,k,n) * EFz; // CGS->MKS convertion of zk
      }
   });
}
#endif
#endif
