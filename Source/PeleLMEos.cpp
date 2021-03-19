#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

void PeleLM::setThermoPress(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::setThermoPress()", setThermoPress);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   for (int lev = 0; lev <= finest_level; ++lev) {
      setThermoPress(lev, a_time);
   }
}

void PeleLM::setThermoPress(int lev, const TimeStamp &a_time) {

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

void PeleLM::calcDivU(int is_init,
                      int computeDiff,
                      int do_avgDown,
                      const TimeStamp &a_time,
                      std::unique_ptr<AdvanceDiffData> &diffData)
{
   BL_PROFILE_VAR("PeleLM::calcDivU()", calcDivU);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   // If requested, compute diffusion terms
   // otherwise assumes it has already been computed and stored in the proper container of diffData
   if (computeDiff) {
      calcDiffusivity(a_time);
      computeDifferentialDiffusionTerms(a_time, diffData, is_init);
   }

   // Compute reactions terms from top level to bottom
   // TODO
   if (is_init) {
   }

   // Assemble divU on each level
   for (int lev = 0; lev <= finest_level; lev++ ) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->divu, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rhoY     = ldata_p->species.const_array(mfi);
         auto const& T        = ldata_p->temp.const_array(mfi);
         auto const& SpecD    = ( a_time == AmrOldTime ) ? diffData->Dn[lev].const_array(mfi,0)
                                                         : diffData->Dnp1[lev].const_array(mfi,0);
         auto const& Fourier  = ( a_time == AmrOldTime ) ? diffData->Dn[lev].const_array(mfi,NUM_SPECIES)
                                                         : diffData->Dnp1[lev].const_array(mfi,NUM_SPECIES);
         auto const& DiffDiff = ( a_time == AmrOldTime ) ? diffData->Dn[lev].const_array(mfi,NUM_SPECIES+1)
                                                         : diffData->Dnp1[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& divu    = ldata_p->divu.array(mfi);
         // TODO reaction
         amrex::ParallelFor(bx, [ rhoY, T, SpecD, Fourier, DiffDiff, divu]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            compute_divu( i, j, k, rhoY, T, SpecD, Fourier, DiffDiff, divu );
         });
      }
   }

   // Average down divU
   if ( do_avgDown ) {
      for (int lev = finest_level; lev > 0; --lev) {
         auto ldataFine_p = getLevelDataPtr(lev,a_time);
         auto ldataCrse_p = getLevelDataPtr(lev-1,a_time);
#ifdef AMREX_USE_EB
         EB_average_down(ldataFine_p->divu,
                         ldataCrse_p->divu,
                         0,1,refRatio(lev-1));
#else
         average_down(ldataFine_p->divu,
                      ldataCrse_p->divu,
                      0,1,refRatio(lev-1));
#endif
      }
   }

   // fillPatch a_time divu to get properly filled ghost cells
   int nGrowDivu = 1;   // TODO: need to make sure it's consistent across
   for (int lev = 0; lev <= finest_level; ++lev) {
      Real time = getTime(lev,a_time);
      auto ldata_p = getLevelDataPtr(lev,a_time);
      fillpatch_divu(lev,time,ldata_p->divu,nGrowDivu);
   }
}

void PeleLM::setTemperature(const TimeStamp &a_time) {
   BL_PROFILE_VAR("PeleLM::setTemperature()", setTemperature);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   for (int lev = 0; lev <= finest_level; ++lev) {
      setTemperature(lev, a_time);
   }
}

void PeleLM::setTemperature(int lev, const TimeStamp &a_time) {

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   auto ldata_p = getLevelDataPtr(lev,a_time);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   {
      for (MFIter mfi(ldata_p->temp,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = ldata_p->density.const_array(mfi);
         auto const& rhoY    = ldata_p->species.const_array(mfi);
         auto const& rhoh    = ldata_p->rhoh.const_array(mfi);
         auto const& T       = ldata_p->temp.array(mfi);

         amrex::ParallelFor(bx, [rho, rhoY, rhoh, T]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getTfromHY( i, j, k, rho, rhoY, rhoh, T);
         });
      }
   }
}

void PeleLM::calc_dPdt(const TimeStamp &a_time,
                       const Vector<MultiFab*> &a_dPdt)
{
   BL_PROFILE_VAR("PeleLM::calc_dPdt()", calc_dPdt);

   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

   for (int lev = 0; lev <= finest_level; ++lev) {
      calc_dPdt(lev, a_time, a_dPdt[lev]);
   }

   // TODO: subcycling version of PeleLM do redistribution when EB

   // Fill ghost cell(s)
   if (a_dPdt[0]->nGrow() > 0) {
      fillpatch_forces(m_cur_time,a_dPdt,a_dPdt[0]->nGrow());
   }
}

void PeleLM::calc_dPdt(int lev,
                       const TimeStamp &a_time,
                       MultiFab* a_dPdt)
{

   auto ldata_p = getLevelDataPtr(lev,a_time);

   Real p_amb = m_pOld;

   if (m_closed_chamber) {
      // TODO
   }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(*a_dPdt,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto const& dPdt  = a_dPdt->array(mfi);
      auto const& P     = ldata_p->rhoRT.const_array(mfi);
      amrex::ParallelFor(bx, [dPdt, P, p_amb, dt=m_dt, dpdt_fac=m_dpdtFactor]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
         dPdt(i,j,k) = (P(i,j,k) - p_amb) / ( dt * P(i,j,k) ) * dpdt_fac;
      });
   }
}
