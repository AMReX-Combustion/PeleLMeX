#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;


void PeleLM::fluxDivergence(Vector<MultiFab> &a_divergence,
                            int div_comp,
                            const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                            int flux_comp,
                            int ncomp,
                            Real scale) {

   BL_PROFILE_VAR("PeleLM::fluxDivergence()", fluxDivergence);
   for (int lev = 0; lev <= finest_level; ++lev) {
      fluxDivergenceLevel(lev,a_divergence[lev], div_comp, a_fluxes[lev], flux_comp, ncomp, scale);
   }
}

void PeleLM::fluxDivergenceLevel(int lev,
                                 MultiFab &a_divergence,
                                 int div_comp,
                                 const Array<MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                                 int flux_comp,
                                 int ncomp,
                                 Real scale) {

   AMREX_ASSERT(a_divergence.nComp() >= div_comp+ncomp);

   // Get the volume
   // TODO: might want to store that somewhere ...
   MultiFab volume(grids[lev], dmap[lev], 1, 0);
   geom[lev].GetVolume(volume);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(a_divergence,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      AMREX_D_TERM(auto const& fluxX = a_fluxes[0]->const_array(mfi,flux_comp);,
                   auto const& fluxY = a_fluxes[1]->const_array(mfi,flux_comp);,
                   auto const& fluxZ = a_fluxes[2]->const_array(mfi,flux_comp););
      auto const& divergence   = a_divergence.array(mfi,div_comp);
      auto const& vol          = volume.const_array(mfi);

#ifdef AMREX_USE_EB
      auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
      auto const& flagfab = ebfactory.getMultiEBCellFlagFab()[mfi];
      auto const& flag    = flagfab.const_array();
#endif

#ifdef AMREX_USE_EB
      if (flagfab.getType(bx) == FabType::covered) {              // Covered boxes
         amrex::ParallelFor(bx, ncomp, [divergence]
         AMREX_GPU_DEVICE( int i, int j, int k, int n) noexcept
         {
            divergence(i,j,k,n) = 0.0;
         });
      } else if (flagfab.getType(bx) != FabType::regular ) {     // EB containing boxes 
         auto vfrac = ebfactory.getVolFrac().const_array(mfi);
         amrex::ParallelFor(bx, [ncomp, flag, vfrac, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            if ( flag(i,j,k).isCovered() ) {
               for (int n = 0; n < nComp; n++) {
                  divergence(i,j,k,n) = 0.0;
               }
            } else if ( flag(i,j,k).isRegular() ) {
               fluxDivergence_K( i, j, k, ncomp,
                                 AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                 vol, scale, divergence);
            } else {
               Real vfracinv = 1.0/vfrac(i,j,k);
               fluxDivergence_K( i, j, k, ncomp,
                                 AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                 vol, scale, divergence);
               for (int n = 0; n < nComp; n++) {
                  divergence(i,j,k,n) *= vfracinv;
               }
            }
         });
      } else {                                                   // Regular boxes
#endif

         amrex::ParallelFor(bx, [ncomp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            fluxDivergence_K( i, j, k, ncomp,
                              AMREX_D_DECL(fluxX, fluxY, fluxZ),
                              vol, scale, divergence);
         });
#ifdef AMREX_USE_EB
      }
#endif
   }

#ifdef AMREX_USE_EB
   {
      MultiFab div_SrcGhostCell(grids,dmap,nComp,a_divergence.nGrow()+2,MFInfo(),Factory());
      div_SrcGhostCell.setVal(0.);
      div_SrcGhostCell.copy(a_divergence, div_comp, 0, ncomp);
      amrex::single_level_weighted_redistribute( {div_SrcGhostCell}, {a_divergence}, *volfrac, div_comp, ncomp, {geom} );
   }
   EB_set_covered(fdiv,0.);
#endif

}
