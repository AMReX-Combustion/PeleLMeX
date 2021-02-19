#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

void PeleLM::setThermoPress(int lev, TimeStamp a_time) {
   BL_PROFILE_VAR("PeleLM::setThermoPress()", setThermoPress);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);
 
   auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(ldata_p->rhoRT,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {    
         const Box& bx = mfi.tilebox();
         auto const& rho     = ldata_p->density.const_array(mfi);
         auto const& rhoY    = ldata_p->species.const_array(mfi);
         auto const& T       = ldata_p->temp.const_array(mfi);
         auto const& P       = ldata_p->rhoRT.array(mfi);

         amrex::ParallelFor(bx, [rho, rhoY, T, P]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {    
            getPGivenRTY( i, j, k, rho, rhoY, T, P ); 
         });  
      }    
   }
}

void PeleLM::calcDivU(int is_init, TimeStamp a_time) {
   BL_PROFILE_VAR("PeleLM::calcDivU()", calcDivU);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   // During initialization, compute transport coefficients
   if ( is_init ) calcDiffusivity(a_time);

   // Compute viscous terms
   int nGrow = 0;                   // No need for ghost cells here
   Vector<MultiFab> viscTerm(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      viscTerm[lev].define(grids[lev], dmap[lev], NUM_SPECIES+2, nGrow);
   }
   computeDifferentialDiffusionTerms(a_time,viscTerm);

   // Compute reactions terms from top level to bottom
   // TODO

   // Assemble divU on each level
   for (int lev = 0; lev <= finest_level; lev++ ) {

      auto ldata_p = getLevelDataPtr(lev,a_time);
      
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->divu, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rhoY    = ldata_p->species.const_array(mfi);
         auto const& T       = ldata_p->temp.const_array(mfi);
         auto const& SpecD   = viscTerm[lev].const_array(mfi,0);
         auto const& Fourier = viscTerm[lev].const_array(mfi,NUM_SPECIES);
         auto const& DiffD   = viscTerm[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& divu    = ldata_p->divu.array(mfi);
         amrex::ParallelFor(bx, [ rhoY, T, SpecD, Fourier, DiffD, divu]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            compute_divu( i, j, k, rhoY, T, SpecD, Fourier, DiffD, divu );
         });
      }
   }
}
