#include <PeleLM.H>
#include <PeleLM_K.H>
#include <reactor.h>

using namespace amrex;

void PeleLM::advanceChemistry(std::unique_ptr<AdvanceAdvData> &advData)
{
   BL_PROFILE_VAR("PeleLM::advanceChemistry()", advanceChemistry);

   for (int lev = 0; lev <= finest_level; ++lev) {
      advanceChemistry(lev, m_dt, advData->Forcing[lev]);
   }
}

void PeleLM::advanceChemistry(int lev,
                              const Real &a_dt,
                              MultiFab &a_extForcing)
{
   auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
   auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
   auto ldataR_p   = getLevelDataReactPtr(lev);

   // TODO Setup covered cells mask
   FabArray<BaseFab<int>> mask(grids[lev],dmap[lev],1,0);
   mask.setVal(1);

   // TODO: try tricks to not do reaction on covered cells

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldataNew_p->density,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx          = mfi.tilebox();
      auto const& rhoY_o     = ldataOld_p->species.const_array(mfi);
      auto const& rhoH_o     = ldataOld_p->rhoh.const_array(mfi);
      auto const& temp_o     = ldataOld_p->temp.const_array(mfi);
      auto const& rhoY_n     = ldataNew_p->species.array(mfi);
      auto const& rhoH_n     = ldataNew_p->rhoh.array(mfi);
      auto const& temp_n     = ldataNew_p->temp.array(mfi);
      auto const& extF_rhoY  = a_extForcing.array(mfi,0);
      auto const& extF_rhoH  = a_extForcing.array(mfi,NUM_SPECIES);
      auto const& fcl        = ldataR_p->functC.array(mfi);
      auto const& mask_arr   = mask.array(mfi);

      // Reset new to old and convert MKS -> CGS
      ParallelFor(bx, [rhoY_o, rhoH_o, temp_o, rhoY_n, rhoH_n, temp_n, extF_rhoY, extF_rhoH]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         for (int n = 0; n < NUM_SPECIES; n++) {
            rhoY_n(i,j,k,n) = rhoY_o(i,j,k,n) * 1.0e-3;
            extF_rhoY(i,j,k,n) *= 1.0e-3;
         }
         temp_n(i,j,k) = temp_o(i,j,k);
         rhoH_n(i,j,k) = rhoH_o(i,j,k) * 10.0;
         extF_rhoH(i,j,k) *= 10.0;
      });

#ifdef AMREX_USE_GPU
      int ncells           = bx.numPts();
      const auto ec = Gpu::ExecutionConfig(ncells);
#endif

      Real dt_incr     = a_dt;
      Real time_chem   = 0;
#ifndef AMREX_USE_GPU
      /* Solve */
      int tmp_fctCn = 0;
      tmp_fctCn = react(bx, rhoY_n, extF_rhoY, temp_n, rhoH_n, extF_rhoH, fcl, mask_arr,
                        dt_incr, time_chem);
      dt_incr   = a_dt;
      time_chem = 0;
#else
      int reactor_type = 2;
      int tmp_fctCn = 0;
      tmp_fctCn = react(bx, rhoY_n, extF_rhoY, temp_n, rhon, extF_rhoH, fcl, mask_arr, 
                        dt_incr, time_chem, reactor_type, amrex::Gpu::gpuStream());
      dt_incr = a_dt;
      time_chem = 0;
#endif

      // Convert CGS -> MKS
      ParallelFor(bx, [rhoY_n, rhoH_n, extF_rhoY, extF_rhoH]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         for (int n = 0; n < NUM_SPECIES; n++) {
            rhoY_n(i,j,k,n) *= 1.0e3;
            extF_rhoY(i,j,k,n) *= 1.0e3;
         }
         rhoH_n(i,j,k) *= 0.1;
         extF_rhoH(i,j,k) *= 0.1;
      });

#ifdef AMREX_USE_GPU
      Gpu::Device::streamSynchronize();
#endif

   }

   // Set reaction term
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldataNew_p->density,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx          = mfi.tilebox();
      auto const& rhoY_o     = ldataOld_p->species.const_array(mfi);
      auto const& rhoY_n     = ldataNew_p->species.const_array(mfi);
      auto const& extF_rhoY  = a_extForcing.const_array(mfi,0);
      auto const& rhoYdot    = ldataR_p->I_R.array(mfi,0);
      Real dt_inv = 1.0/a_dt;
      ParallelFor(bx, NUM_SPECIES, [rhoY_o, rhoY_n, extF_rhoY, rhoYdot, dt_inv]
      AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
      {
         rhoYdot(i,j,k,n) = - ( rhoY_o(i,j,k,n) - rhoY_n(i,j,k,n) ) * dt_inv - extF_rhoY(i,j,k,n);
      });
   }
}

void PeleLM::computeInstantaneousReactionRate(const Vector<MultiFab*> &I_R,
                                              const TimeStamp &a_time)
{
   BL_PROFILE_VAR("PeleLM::computeInstantaneousReactionRate()", computeInstantaneousReactionRate);

   for (int lev = 0; lev <= finest_level; ++lev) {
      // TODO Setup covered cells mask
      MultiFab mask(grids[lev],dmap[lev],1,0);
      mask.setVal(1.0);
      computeInstantaneousReactionRate(lev, a_time, mask, I_R[lev]);
   }
}

void PeleLM::computeInstantaneousReactionRate(int lev,
                                              const TimeStamp &a_time,
                                              const MultiFab &a_mask,
                                              MultiFab* a_I_R)
{
   auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldata_p->species,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& rhoY    = ldata_p->species.const_array(mfi);
      auto const& rhoH    = ldata_p->rhoh.const_array(mfi);
      auto const& T       = ldata_p->temp.const_array(mfi);
      auto const& mask    = a_mask.const_array(mfi);
      auto const& rhoYdot = a_I_R->array(mfi);

      amrex::ParallelFor(bx, [rhoY, rhoH, T, mask, rhoYdot]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         reactionRateRhoY( i, j, k, rhoY, rhoH, T, mask,
                           rhoYdot );
      });
   }
}

void PeleLM::getScalarReactForce(std::unique_ptr<AdvanceAdvData> &advData)
{
   // The differentialDiffusionUpdate just provided the {np1,kp1} AD state
   // -> use it to build the external forcing for the chemistry
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get t^{n} t^{np1} data pointer
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
      auto ldataR_p = getLevelDataReactPtr(lev);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(advData->Forcing[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rhoY_o     = ldataOld_p->species.const_array(mfi);
         auto const& rhoH_o     = ldataOld_p->rhoh.const_array(mfi);
         auto const& rhoY_n     = ldataNew_p->species.const_array(mfi);
         auto const& rhoH_n     = ldataNew_p->rhoh.const_array(mfi);
         auto const& react      = ldataR_p->I_R.const_array(mfi,0);
         auto const& extF_rhoY  = advData->Forcing[lev].array(mfi,0);
         auto const& extF_rhoH  = advData->Forcing[lev].array(mfi,NUM_SPECIES);
         amrex::Real dtinv   = 1.0/m_dt;
         amrex::ParallelFor(bx, [rhoY_o, rhoH_o, rhoY_n, rhoH_n, react,
                                 extF_rhoY, extF_rhoH, dtinv]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            for (int n = 0; n < NUM_SPECIES; n++) {
               extF_rhoY(i,j,k,n) = (rhoY_n(i,j,k,n) - rhoY_o(i,j,k,n)) * dtinv - react(i,j,k,n);
            }
            extF_rhoH(i,j,k) = (rhoH_n(i,j,k) - rhoH_o(i,j,k)) * dtinv;
         });
      }
   }
}
