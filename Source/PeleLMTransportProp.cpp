#include <PeleLM.H>
#include <PeleLM_K.H>
#include <pelelm_prob.H>
#ifdef PELE_USE_EFIELD
#include <PeleLMEF_K.H>
#endif

using namespace amrex;

void PeleLM::calcViscosity(const TimeStamp &a_time) {
   BL_PROFILE("PeleLM::calcViscosity()");

   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

      if (m_incompressible) {
         ldata_p->visc_cc.setVal(m_mu);
      } else {

         // Transport data pointer
         auto const* ltransparm = trans_parms.device_trans_parm();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldata_p->visc_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& gbx     = mfi.growntilebox();
            auto const& rhoY   = ldata_p->state.const_array(mfi,FIRSTSPEC);
            auto const& T      = ldata_p->state.const_array(mfi,TEMP);
            auto const& mu     = ldata_p->visc_cc.array(mfi,0);

            amrex::ParallelFor(gbx, [rhoY, T, mu, ltransparm]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               getVelViscosity( i, j, k, rhoY, T, mu, ltransparm);
            });
         }
      }
   }
}

void PeleLM::calcDiffusivity(const TimeStamp &a_time) {
   BL_PROFILE("PeleLM::calcDiffusivity()");

   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,a_time);

      // Transport data pointer
      auto const* ltransparm = trans_parms.device_trans_parm();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->diff_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx     = mfi.growntilebox();
         auto const& rhoY   = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& T      = ldata_p->state.const_array(mfi,TEMP);
         auto const& rhoD   = ldata_p->diff_cc.array(mfi,0);
         auto const& lambda = ldata_p->diff_cc.array(mfi,NUM_SPECIES);
         auto const& mu     = ldata_p->diff_cc.array(mfi,NUM_SPECIES+1);

         // TODO: unity Lewis

         amrex::ParallelFor(gbx, [rhoY, T, rhoD, lambda, mu, ltransparm]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getTransportCoeff( i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
         });
#ifdef PELE_USE_EFIELD
         auto const& Ks   = ldata_p->mob_cc.array(mfi,0);
         auto eos = pele::physics::PhysicsType::eos();
         Real mwt[NUM_SPECIES] = {0.0};
         eos.molecular_weight(mwt);
         amrex::ParallelFor(gbx, [rhoY, rhoD, T, Ks, mwt, zk=zk]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getKappaSp( i, j, k, mwt, zk, rhoY, rhoD, T, Ks);
         });
#endif
      }

   }
}

Array<MultiFab,AMREX_SPACEDIM>
PeleLM::getDiffusivity(int lev, int beta_comp, int ncomp, int doZeroVisc,
                       Vector<BCRec> bcrec,
                       MultiFab const& beta_cc)
{
   BL_PROFILE("PeleLM::getDiffusivity()");

   AMREX_ASSERT(bcrec.size() >= ncomp);
   AMREX_ASSERT(beta_cc.nComp() >= beta_comp+ncomp);

   const auto& ba = beta_cc.boxArray();
   const auto& dm = beta_cc.DistributionMap();
   const auto& factory = beta_cc.Factory();
   Array<MultiFab,AMREX_SPACEDIM> beta_ec{AMREX_D_DECL(MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                                                dm, ncomp, 0, MFInfo(), factory),
                                                       MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                                                dm, ncomp, 0, MFInfo(), factory),
                                                       MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                                                dm, ncomp, 0, MFInfo(), factory))};

   const Box& domain = geom[lev].Domain();

#ifdef AMREX_USE_EB
   // EB : use EB CCentroid -> FCentroid
   EB_interp_CellCentroid_to_FaceCentroid(beta_cc, GetArrOfPtrs(beta_ec), beta_comp, 0, ncomp, geom[lev], bcrec);
   EB_set_covered_faces(GetArrOfPtrs(beta_ec),1.234e40);
#else
   // NON-EB : use cen2edg_cpp
   bool use_harmonic_avg = m_harm_avg_cen2edge ? true : false;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(beta_cc,TilingIfNotGPU()); mfi.isValid();++mfi)
   {
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++)
      {
         const Box ebx = mfi.nodaltilebox(idim);
         const Box& edomain = amrex::surroundingNodes(domain,idim);
         const auto& diff_c  = beta_cc.const_array(mfi,beta_comp);
         const auto& diff_ec = beta_ec[idim].array(mfi);
         const auto bc_lo = bcrec[0].lo(idim);
         const auto bc_hi = bcrec[0].hi(idim);
         amrex::ParallelFor(ebx, [idim, ncomp, bc_lo, bc_hi, use_harmonic_avg, diff_c, diff_ec, edomain]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            int idx[3] = {i,j,k};
            bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) && ( idx[idim] <= edomain.smallEnd(idim) ) );
            bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) && ( idx[idim] >= edomain.bigEnd(idim) ) );
            cen2edg_cpp( i, j, k, idim, ncomp, use_harmonic_avg, on_lo, on_hi, diff_c, diff_ec);
         });
      }
   }
#endif

   // Enable zeroing diffusivity on faces to produce walls
   if (doZeroVisc) {
      const auto geomdata = geom[lev].data();
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         const Box& edomain = amrex::surroundingNodes(domain,idim);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(beta_ec[idim],TilingIfNotGPU()); mfi.isValid();++mfi) {
            const Box ebx = mfi.tilebox();
            const auto& diff_ec = beta_ec[idim].array(mfi);
            amrex::ParallelFor(ebx, [=]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                zero_visc(i, j, k, diff_ec, geomdata, edomain, idim, beta_comp, ncomp);
            });
         }
      }
   }

   return beta_ec;
}
