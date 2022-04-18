#include <PeleLM.H>
#include <PeleLM_K.H>
#include <hydro_utils.H>

using namespace amrex;

void PeleLM::fluxDivergence(const Vector<MultiFab*> &a_divergence,
                            int div_comp,
                            const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                            int flux_comp,
                            int ncomp,
                            int intensiveFluxes,
                            Real scale) {

   BL_PROFILE("PeleLM::fluxDivergence()");
   if (intensiveFluxes) {        // Fluxes are intensive -> need area scaling in div
      for (int lev = 0; lev <= finest_level; ++lev) {
         intFluxDivergenceLevel(lev,*a_divergence[lev], div_comp, a_fluxes[lev], flux_comp,
                                ncomp, scale);
      }
   } else {                      // Fluxes are extensive
      for (int lev = 0; lev <= finest_level; ++lev) {
         extFluxDivergenceLevel(lev,*a_divergence[lev], div_comp, a_fluxes[lev], flux_comp,
                                ncomp, scale);
      }
   }
}

void PeleLM::fluxDivergenceRD(const Vector<const MultiFab*> &a_state,
                              int state_comp,
                              const Vector<MultiFab*> &a_divergence,
                              int div_comp,
                              const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                              int flux_comp,
                              int ncomp,
                              int intensiveFluxes,
                              const BCRec *state_bc_d,
                              const Real &scale,
                              const Real &a_dt)
{
    BL_PROFILE("PeleLM::fluxDivergenceRD()");
#ifdef AMREX_USE_EB
    for (int lev = 0; lev <= finest_level; ++lev) {
        //----------------------------------------------------------------
        // Use a temporary MF to hold divergence before redistribution
        int nGrow_divTmp= 3;
        MultiFab divTmp(grids[lev],dmap[lev],ncomp,nGrow_divTmp,MFInfo(),EBFactory(lev));
        divTmp.setVal(0.0);
        if (intensiveFluxes) {        // Fluxes are intensive -> need area scaling in div
            intFluxDivergenceLevel(lev, divTmp, 0, a_fluxes[lev], flux_comp, ncomp, scale);
        } else {                      // Fluxes are extensive
            extFluxDivergenceLevel(lev, divTmp, 0, a_fluxes[lev], flux_comp, ncomp, scale);
        }
        
        // Need FillBoundary before redistribution
        divTmp.FillBoundary(geom[lev].periodicity());
        
        // Redistribute diffusion term
        redistributeDiff(lev, a_dt,
                         divTmp, 0,
                         *a_divergence[lev], div_comp,
                         *a_state[lev], state_comp,
                         ncomp,
                         state_bc_d,
                         geom[lev]);
    }
#else
    amrex::ignore_unused(a_state);
    amrex::ignore_unused(state_comp);
    amrex::ignore_unused(state_bc_d);
    amrex::ignore_unused(a_dt);
    fluxDivergence(a_divergence, div_comp, a_fluxes, flux_comp, ncomp, intensiveFluxes, scale);
#endif
}

void PeleLM::extFluxDivergenceLevel(int lev,
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

#ifdef AMREX_USE_EB
   auto const& ebfact = EBFactory(lev);
#endif

#ifdef AMREX_USE_OMP
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
      auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
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
         auto vfrac = ebfact.getVolFrac().const_array(mfi);
         amrex::ParallelFor(bx, [ncomp, flag, vfrac, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            if ( flag(i,j,k).isCovered() ) {
               for (int n = 0; n < ncomp; n++) {
                  divergence(i,j,k,n) = 0.0;
               }
            } else if ( flag(i,j,k).isRegular() ) {
               extFluxDivergence_K( i, j, k, ncomp,
                                    AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                    vol, scale, divergence);
            } else {
               Real vfracinv = 1.0/vfrac(i,j,k);
               extFluxDivergence_K( i, j, k, ncomp,
                                    AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                    vol, scale, divergence);
               for (int n = 0; n < ncomp; n++) {
                  divergence(i,j,k,n) *= vfracinv;
               }
            }
         });
      } else {                                                   // Regular boxes
#endif

         amrex::ParallelFor(bx, [ncomp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            extFluxDivergence_K( i, j, k, ncomp,
                                 AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                 vol, scale, divergence);
         });
#ifdef AMREX_USE_EB
      }
#endif
   }
}

void PeleLM::intFluxDivergenceLevel(int lev,
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

   // Get area
   const Real* dx = Geom(lev).CellSize();
#if ( AMREX_SPACEDIM == 2 )
   Real areax = dx[1];
   Real areay = dx[0];
#elif ( AMREX_SPACEDIM == 3 )
   Real areax = dx[1]*dx[2];
   Real areay = dx[0]*dx[2];
   Real areaz = dx[0]*dx[1];
#endif

   // Get areafrac if EB
#ifdef AMREX_USE_EB
      auto const& ebfact = EBFactory(lev);
      Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
      areafrac  = ebfact.getAreaFrac();
#endif

#ifdef AMREX_USE_OMP
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
      auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
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
         auto vfrac = ebfact.getVolFrac().const_array(mfi);
         AMREX_D_TERM( const auto& afrac_x = areafrac[0]->array(mfi);,
                       const auto& afrac_y = areafrac[1]->array(mfi);,
                       const auto& afrac_z = areafrac[2]->array(mfi););
         amrex::ParallelFor(bx, [ncomp, flag, vfrac, divergence,
                                 AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                 AMREX_D_DECL(afrac_x, afrac_y, afrac_z),
                                 AMREX_D_DECL(areax, areay, areaz),
                                 vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            if ( flag(i,j,k).isCovered() ) {
               for (int n = 0; n < ncomp; n++) {
                  divergence(i,j,k,n) = 0.0;
               }
            } else if ( flag(i,j,k).isRegular() ) {
               intFluxDivergence_K( i, j, k, ncomp,
                                    AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                    AMREX_D_DECL(areax, areay, areaz),
                                    vol, scale, divergence);
            } else {
               Real vfracinv = 1.0/vfrac(i,j,k);
               EB_intFluxDivergence_K( i, j, k, ncomp,
                                    AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                    AMREX_D_DECL(afrac_x, afrac_y, afrac_z),
                                    AMREX_D_DECL(areax, areay, areaz),
                                    vol, scale, divergence);
               for (int n = 0; n < ncomp; n++) {
                  divergence(i,j,k,n) *= vfracinv;
               }
            }
         });
      } else {                                                   // Regular boxes
#endif

         amrex::ParallelFor(bx, [ncomp, divergence,
                                 AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                 AMREX_D_DECL(areax, areay, areaz),
                                 vol, scale]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            intFluxDivergence_K( i, j, k, ncomp,
                                 AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                 AMREX_D_DECL(areax, areay, areaz),
                                 vol, scale, divergence);
         });
#ifdef AMREX_USE_EB
      }
#endif
   }
}

void PeleLM::advFluxDivergence(int a_lev,
                               MultiFab &a_divergence, int div_comp,
                               MultiFab &a_divu,
                               const Array<const MultiFab*,AMREX_SPACEDIM> &a_fluxes, int flux_comp,
                               const Array<const MultiFab*,AMREX_SPACEDIM> &a_faceState, int face_comp,
                               int ncomp,
                               int const* l_conserv_d,
                               const Geometry &a_geom,
                               amrex::Real scale,
                               bool fluxes_are_area_weighted)
{
    BL_PROFILE("PeleLM::advFluxDivergence()");
    AMREX_ASSERT(a_divergence.nComp() >= div_comp+ncomp);
    AMREX_ASSERT(a_fluxes[0]->nComp() >= flux_comp+ncomp);
    AMREX_ASSERT(a_faceState[0]->nComp() >= face_comp+ncomp);

#ifdef AMREX_USE_EB
    auto const& ebfact = EBFactory(a_lev);
#else
    amrex::ignore_unused(a_lev); 
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_divergence,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    
        Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
#endif
    
        // Get the divergence
        auto const& div_arr  = a_divergence.array(mfi,div_comp);
        AMREX_D_TERM(auto const& fx = a_fluxes[0]->const_array(mfi,flux_comp);,
                     auto const& fy = a_fluxes[1]->const_array(mfi,flux_comp);,
                     auto const& fz = a_fluxes[2]->const_array(mfi,flux_comp);)

#ifdef AMREX_USE_EB
        auto const& vfrac_arr = ebfact.getVolFrac().const_array(mfi);
        if (flagfab.getType(bx) != FabType::covered) {
           HydroUtils::EB_ComputeDivergence(bx, div_arr,
                                            AMREX_D_DECL(fx,fy,fz),
                                            vfrac_arr,
                                            ncomp, a_geom,
                                            scale, fluxes_are_area_weighted);
        }
#else
        HydroUtils::ComputeDivergence(bx, div_arr,
                                      AMREX_D_DECL(fx,fy,fz),
                                      ncomp, a_geom,
                                      scale, fluxes_are_area_weighted);
#endif

        // If convective, we define u dot grad q = div (u q) - q div(u)
        // averaging face and t^{n+1/2} q to the cell center
        auto const& divu_arr    = a_divu.const_array(mfi);
        AMREX_D_TERM(auto const& facex = a_faceState[0]->const_array(mfi,face_comp);,
                     auto const& facey = a_faceState[1]->const_array(mfi,face_comp);,
                     auto const& facez = a_faceState[2]->const_array(mfi,face_comp);)

        bool regular = true;
#ifdef AMREX_USE_EB
        regular = (flagfab.getType(bx) == FabType::regular);
#endif
        if (regular) {
            ParallelFor(bx, ncomp, [=]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
                if (!l_conserv_d[n]) {
                    Real qavg  = AMREX_D_TERM(  facex(i,j,k,n) + facex(i+1,j,k,n),
                                              + facey(i,j,k,n) + facey(i,j+1,k,n),
                                              + facez(i,j,k,n) + facez(i,j,k+1,n));
#if (AMREX_SPACEDIM == 2)
                    qavg *= 0.25;
#else
                    qavg /= 6.0;
#endif
                    // Note that because we define adv update as MINUS div(u q), here we add q div (u)
                    div_arr(i,j,k,n) += qavg*divu_arr(i,j,k);
                }
            });
        }
#ifdef AMREX_USE_EB
        else {
            if (flagfab.getType(bx) == FabType::covered) {
                AMREX_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n, {div_arr(i,j,k,n) = 0.0;});
            } else {
                AMREX_D_TERM(auto const& apx_arr = ebfact.getAreaFrac()[0]->const_array(mfi);,
                             auto const& apy_arr = ebfact.getAreaFrac()[1]->const_array(mfi);,
                             auto const& apz_arr = ebfact.getAreaFrac()[2]->const_array(mfi););
                ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    if (!l_conserv_d[n] && vfrac_arr(i,j,k) > 0.) {
                        Real qwsum  = AMREX_D_TERM(  apx_arr(i,j,k) * facex(i,j,k,n) + apx_arr(i+1,j,k) * facex(i+1,j,k,n),
                                                   + apy_arr(i,j,k) * facey(i,j,k,n) + apy_arr(i,j+1,k) * facey(i,j+1,k,n),
                                                   + apz_arr(i,j,k) * facez(i,j,k,n) + apz_arr(i,j,k+1) * facez(i,j,k+1,n));
                        Real areasum  = AMREX_D_TERM(  apx_arr(i,j,k) + apx_arr(i+1,j,k),
                                                     + apy_arr(i,j,k) + apy_arr(i,j+1,k),
                                                     + apz_arr(i,j,k) + apz_arr(i,j,k+1));
                        // Note that because we define adv update as MINUS div(u q), here we add q div (u)
                        div_arr(i,j,k,n) += qwsum / areasum * divu_arr(i,j,k);
                    }
                });
            }
        }
#endif
    }
}

void
PeleLM::floorSpecies(const TimeStamp &a_time)
{
   AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);
   if (!m_floor_species) {
      return;
   }
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,a_time);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rhoY    = ldata_p->state.array(mfi,FIRSTSPEC);
         amrex::ParallelFor(bx, [rhoY]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            fabMinMax( i, j, k, NUM_SPECIES, 0.0, AMREX_REAL_MAX, rhoY);
         });
#ifdef PELE_USE_EFIELD
         auto const& nE    = ldata_p->nE.array(mfi);
         amrex::ParallelFor(bx, [nE]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            fabMinMax( i, j, k, 1, 0.0, AMREX_REAL_MAX, nE);
         });
#endif
      }
   }
}

void PeleLM::resetCoveredMask()
{
   if (m_resetCoveredMask) {

      if (m_verbose) Print() << " Resetting fine-covered cells mask \n";

      for (int lev = 0; lev < finest_level; ++lev) {
         BoxArray baf = grids[lev+1];
         baf.coarsen(ref_ratio[lev]);
         m_coveredMask[lev]->setVal(1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         {
            std::vector< std::pair<int,Box> > isects;
            for (MFIter mfi(*m_coveredMask[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               auto const& mask = m_coveredMask[lev]->array(mfi);
               baf.intersections(grids[lev][mfi.index()],isects);
               for (const auto& is : isects)
               {
                  amrex::ParallelFor(is.second, [mask]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     mask(i,j,k) = 0;
                  });
               }
            }
         }

         //----------------------------------------------------------------------------
         // Setup a BoxArray for the chemistry

         // Get an uncovered BoxArray
         BoxArray baCompDom = complementIn(geom[lev].Domain(), baf);
         BoxArray baUnCovered = intersect(baCompDom,grids[lev]);

         // Assemble a BoxArray with covered and uncovered ones + flags
         BoxList bl(grids[lev].ixType());
         bl.reserve(baUnCovered.size()+baf.size());
         m_baChemFlag[lev].resize(baUnCovered.size()+baf.size());
         int bxIdx = 0;
         for (int i = 0, n = baUnCovered.size(); i < n; ++i) {  // Append uncovered boxes
            bl.push_back(baUnCovered[i]);
            m_baChemFlag[lev][bxIdx] = 1;
            bxIdx += 1;
         }
         for (int i = 0, n = baf.size(); i < n; ++i) {          // Append covered boxes
            bl.push_back(baf[i]);
            m_baChemFlag[lev][bxIdx] = 0;
            bxIdx += 1;
         }
         m_baChem[lev].reset(new BoxArray(std::move(bl)));
         m_dmapChem[lev].reset(new DistributionMapping(*m_baChem[lev]));
      }

      // Switch off trigger
      m_resetCoveredMask = 0;
   }

   //----------------------------------------------------------------------------
   // Need to compute the uncovered volume
   // TODO Might need to recompute if new levels are added with EB.
   if (m_uncoveredVol < 0.0 ) {
      Vector<MultiFab> dummy(finest_level+1);
      for (int lev = 0; lev <= finest_level; ++lev) {
         dummy[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *m_factory[lev]);
         dummy[lev].setVal(1.0);
      }
      m_uncoveredVol = MFSum(GetVecOfConstPtrs(dummy),0);
   }

}

// Return a unique_ptr with the entire derive
std::unique_ptr<MultiFab>
PeleLM::derive(const std::string &a_name,
               Real               a_time,
               int                lev,
               int                nGrow)
{
   AMREX_ASSERT(nGrow >= 0);

   std::unique_ptr<MultiFab> mf;

   bool itexists = derive_lst.canDerive(a_name) || isStateVariable(a_name);

   if ( !itexists ) {
      amrex::Error("PeleLM::derive(): unknown variable: "+a_name);
   }

   const PeleLMDeriveRec* rec = derive_lst.get(a_name);

   if (rec) {        // This is a derived variable
      mf.reset(new MultiFab(grids[lev], dmap[lev], rec->numDerive(), nGrow));
      std::unique_ptr<MultiFab> statemf = fillPatchState(lev, a_time, m_nGrowState);
      // Get pressure: TODO no fillpatch for pressure just yet, simply get new state
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      auto stateBCs = fetchBCRecArray(VELX,NVAR);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const Box& bx = mfi.growntilebox(nGrow);
          FArrayBox& derfab = (*mf)[mfi];
          FArrayBox const& statefab = (*statemf)[mfi];
          FArrayBox const& pressfab = ldata_p->press[mfi];
          rec->derFunc()(bx, derfab, 0, rec->numDerive(), statefab, pressfab, geom[lev], a_time, stateBCs, lev);
      }
   } else {          // This is a state variable
      mf.reset(new MultiFab(grids[lev], dmap[lev], 1, nGrow));
      int idx = stateVariableIndex(a_name);
      std::unique_ptr<MultiFab> statemf = fillPatchState(lev, a_time, nGrow);
      MultiFab::Copy(*mf,*statemf,idx,0,1,nGrow);
   }

   return mf;
}

// Return a unique_ptr with only the required component of a derive
std::unique_ptr<MultiFab>
PeleLM::deriveComp(const std::string &a_name,
                   Real               a_time,
                   int                lev,
                   int                nGrow)
{
   AMREX_ASSERT(nGrow >= 0);

   std::unique_ptr<MultiFab> mf;

   bool itexists = derive_lst.canDerive(a_name) || isStateVariable(a_name);

   if ( !itexists ) {
      amrex::Error("PeleLM::derive(): unknown variable: "+a_name);
   }

   const PeleLMDeriveRec* rec = derive_lst.get(a_name);

   if (rec) {        // This is a derived variable
      mf.reset(new MultiFab(grids[lev], dmap[lev], 1, nGrow));
      std::unique_ptr<MultiFab> statemf = fillPatchState(lev, a_time, m_nGrowState);
      // Get pressure: TODO no fillpatch for pressure just yet, simply get new state
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      auto stateBCs = fetchBCRecArray(VELX,NVAR);

      // Temp MF for all the derive components
      MultiFab derTemp(grids[lev], dmap[lev], rec->numDerive(), nGrow);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*mf,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
          const Box& bx = mfi.growntilebox(nGrow);
          FArrayBox& derfab = derTemp[mfi];
          FArrayBox const& statefab = (*statemf)[mfi];
          FArrayBox const& pressfab = ldata_p->press[mfi];
          rec->derFunc()(bx, derfab, 0, rec->numDerive(), statefab, pressfab, geom[lev], a_time, stateBCs, lev);
      }
      // Copy into outgoing unique_ptr
      int derComp = rec->variableComp(a_name);
      if ( derComp < 0 ) {
         amrex::Error("PeleLM::deriveComp(): unknown derive component: " + a_name + " of " + rec->variableName(1000));
      }
      MultiFab::Copy(*mf,derTemp,derComp,0,1,nGrow);
   } else {          // This is a state variable
      mf.reset(new MultiFab(grids[lev], dmap[lev], 1, nGrow));
      int idx = stateVariableIndex(a_name);
      std::unique_ptr<MultiFab> statemf = fillPatchState(lev, a_time, nGrow);
      MultiFab::Copy(*mf,*statemf,idx,0,1,nGrow);
   }

   return mf;
}

Real
PeleLM::MLNorm0(const Vector<const MultiFab*> &a_MF)
{
   Real r = 0.0;
   for (int lev = 0; lev < a_MF.size(); ++lev) {
      if (lev != finest_level) {
         r = std::max(r, a_MF[lev]->norm0(*m_coveredMask[lev],0,0,true));
      } else {
         r = std::max(r, a_MF[lev]->norm0(0,0,true,true));
      }
   }
   ParallelDescriptor::ReduceRealMax(r);
   return r;
}

Vector<Real>
PeleLM::MLNorm0(const Vector<const MultiFab*> &a_MF,
                int startcomp, int ncomp)
{
   AMREX_ASSERT(a_MF[0]->nComp() >= startcomp+ncomp);
   Vector<Real> r(ncomp);
   for (int n = 0; n < ncomp; n++) {
      r[n] = 0.0;
   }
   for (int lev = 0; lev < a_MF.size(); ++lev) {
      if (lev != finest_level) {
         for (int n = 0; n < ncomp; n++) {
            r[n] = std::max(r[n], a_MF[lev]->norm0(*m_coveredMask[lev],startcomp+n,0,true));
         }
      } else {
         for (int n = 0; n < ncomp; n++) {
            r[n] = std::max(r[n], a_MF[lev]->norm0(startcomp+n,0,true,true));
         }
      }
   }
   ParallelDescriptor::ReduceRealMax(r.data(),ncomp);
   return r;
}

bool
PeleLM::isStateVariable(const std::string &a_name)
{
   for (std::list<std::tuple<int,std::string>>::const_iterator li = stateComponents.begin(),
        End = stateComponents.end(); li != End; ++li)
   {
      if (std::get<1>(*li) == a_name) {
         return true;
      }
   }
   return false;
}

int
PeleLM::stateVariableIndex(const std::string &a_name)
{
   int idx = -1;
   if (!isStateVariable(a_name)) {
      amrex::Error("PeleLM::stateVariableIndex(): unknown State variable: "+a_name);
   }
   for (std::list<std::tuple<int,std::string>>::const_iterator li = stateComponents.begin(),
        End = stateComponents.end(); li != End; ++li)
   {
      if (std::get<1>(*li) == a_name) {
         idx = std::get<0>(*li);
      }
   }
   return idx;
}

Vector<int>
PeleLM::fetchAdvTypeArray(int scomp, int ncomp)
{
   Vector<int> types(ncomp);
   for (int comp = 0; comp < ncomp; comp++) {
      types[comp] = m_AdvTypeState[scomp+comp];
   }
   return types;
}

Vector<int>
PeleLM::fetchDiffTypeArray(int scomp, int ncomp)
{
   Vector<int> types(ncomp);
   for (int comp = 0; comp < ncomp; comp++) {
      types[comp] = m_DiffTypeState[scomp+comp];
   }
   return types;
}

Real 
PeleLM::MFSum (const Vector<const MultiFab*> &a_mf, int comp)
{
    // Get the integral of the MF, not including the fine-covered cells

    Real  volwgtsum = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const Real* dx = geom[lev].CellSize();

       // Use amrex::ReduceSum
       Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);

#ifdef AMREX_USE_EB
       auto const& ebfact = dynamic_cast<EBFArrayBoxFactory const&>(Factory(lev));
       auto const& vfrac = ebfact.getVolFrac();
   
       Real sm = 0.0;
       if ( lev != finest_level ) {
          sm = amrex::ReduceSum(*a_mf[lev], vfrac, *m_coveredMask[lev], 0, [vol, comp]
          AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, 
                                                Array4<Real const> const& vf_arr,
                                                Array4<int  const> const& covered_arr) -> Real
          {
              Real sum = 0.0;
              AMREX_LOOP_3D(bx, i, j, k,
              {
                  sum += mf_arr(i,j,k,comp) * vf_arr(i,j,k) * vol * static_cast<Real>(covered_arr(i,j,k));
              });
              return sum;
          });
       } else {
          sm = amrex::ReduceSum(*a_mf[lev], vfrac, 0, [vol, comp]
          AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, 
                                                Array4<Real const> const& vf_arr) -> Real
          {
              Real sum = 0.0;
              AMREX_LOOP_3D(bx, i, j, k,
              {
                  sum += mf_arr(i,j,k,comp) * vf_arr(i,j,k) * vol;
              });
              return sum;
          });
       }
#else
       Real sm = 0.0;
       if ( lev != finest_level ) {
          sm = amrex::ReduceSum(*a_mf[lev], *m_coveredMask[lev], 0, [vol, comp]
          AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr, Array4<int const> const& covered_arr) -> Real
          {
              Real sum = 0.0;
              AMREX_LOOP_3D(bx, i, j, k,
              {
                  sum += mf_arr(i,j,k,comp) * vol * static_cast<Real>(covered_arr(i,j,k));
              });
              return sum;
          });
       } else {
          sm = amrex::ReduceSum(*a_mf[lev], 0, [vol, comp]
          AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& mf_arr) -> Real
          {
              Real sum = 0.0;
              AMREX_LOOP_3D(bx, i, j, k,
              {
                  sum += mf_arr(i,j,k,comp) * vol;
              });
              return sum;
          });
       }
#endif

        volwgtsum += sm;
    } // lev

    ParallelDescriptor::ReduceRealSum(volwgtsum);

    return volwgtsum;
}

/*
Array<Real,3>
PeleLM::MFStat (const Vector<const MultiFab*> &a_mf, int comp)
{
   // Get the min/max/mean of a given component, not including the fine-covered cells
}
*/

void PeleLM::setTypicalValues(const TimeStamp &a_time, int is_init)
{
    // Get state Max/Min
    auto stateMax = ( m_incompressible ) ? MLmax(GetVecOfConstPtrs(getStateVect(a_time)),0,AMREX_SPACEDIM)
                                         : MLmax(GetVecOfConstPtrs(getStateVect(a_time)),0,NVAR);
    auto stateMin = ( m_incompressible ) ? MLmin(GetVecOfConstPtrs(getStateVect(a_time)),0,AMREX_SPACEDIM)
                                         : MLmin(GetVecOfConstPtrs(getStateVect(a_time)),0,NVAR);

    // Fill typical values vector
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
       typical_values[idim] = std::max(stateMax[VELX+idim],std::abs(stateMin[VELX+idim]));
    }
    
    if (!m_incompressible) {
        // Average between max/min
        typical_values[DENSITY] = 0.5 * (stateMax[DENSITY] + stateMin[DENSITY]);
        for (int n = 0; n < NUM_SPECIES; n++) {
            typical_values[FIRSTSPEC+n] = 0.5 * (stateMax[FIRSTSPEC+n] + stateMin[FIRSTSPEC+n]) / typical_values[DENSITY];
        }
        typical_values[RHOH] = 0.5 * (stateMax[RHOH] + stateMin[RHOH]) / typical_values[DENSITY];
        typical_values[TEMP] = 0.5 * (stateMax[TEMP] + stateMin[TEMP]);
        typical_values[RHORT] = m_pOld;

        // Pass into chemsitry if requested
        updateTypicalValuesChem();
    }

    if (is_init || m_verbose > 1) {
        std::string PrettyLine = std::string(78, '=') + "\n";
        Print() << PrettyLine;
        Print() << "Typical values: " << '\n';
        Print() << "\tVelocity: ";
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            Print() << typical_values[idim] << ' ';
        }
        Print() << '\n';
        if (!m_incompressible) {
            Print() << "\tDensity: " << typical_values[DENSITY] << '\n';
            Print() << "\tTemp:    " << typical_values[TEMP]    << '\n';
            Print() << "\tH:       " << typical_values[RHOH]    << '\n';
            Vector<std::string> spec_names;
            pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
            for (int n = 0; n < NUM_SPECIES; n++) {
                Print() << "\tY_" << spec_names[n] << ": " << typical_values[FIRSTSPEC+n] <<'\n';
            }
        }
        Print() << PrettyLine;
    }
}

void PeleLM::updateTypicalValuesChem()
{
    if (m_useTypValChem && m_do_react) {
        if (m_verbose > 2) Print() << " Update chemistry typical values \n";
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            Vector<Real> typical_values_chem;
            typical_values_chem.resize(NUM_SPECIES+1);
            for (int i=0; i<NUM_SPECIES; ++i) {
              typical_values_chem[i] =
                amrex::max(m_typicalYvalMin,
                           typical_values[FIRSTSPEC+i] * typical_values[DENSITY] * 1.E-3); // CGS -> MKS conversion
            }
            typical_values_chem[NUM_SPECIES] = typical_values[TEMP];
            m_reactor->set_typ_vals_ode(typical_values_chem);
        }
    }
}

// MultiFab max, exlucing EB-covered/fine-covered cells, local
Real
PeleLM::MFmax(const MultiFab *a_MF,
              const iMultiFab &a_mask,
              int comp)
{
    Real mx = std::numeric_limits<Real>::lowest();

#ifdef AMREX_USE_EB
    if ( a_MF->hasEBFabFactory() )
    {    
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_MF->Factory());
        auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& flagsma = flags.const_arrays();
            auto const& ma = a_MF->const_arrays();
            auto const& mask = a_mask.const_arrays();
            mx = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, *a_MF, IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {    
                if (flagsma[box_no](i,j,k).isCovered() ||
                    !mask[box_no](i,j,k)) {
                    return AMREX_REAL_LOWEST;
                } else {
                    return ma[box_no](i,j,k,comp);
                }    
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:mx)
#endif
            for (MFIter mfi(*a_MF,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                if (flags[mfi].getType(bx) != FabType::covered) {
                    auto const& flag = flags.const_array(mfi);
                    auto const& a = a_MF->const_array(mfi);
                    auto const& mask = a_mask.const_array(mfi);
                    AMREX_LOOP_3D(bx, i, j, k,
                    {
                        if (!flag(i,j,k).isCovered() && mask(i,j,k)) {
                            mx = std::max(mx, a(i,j,k,comp));
                        }
                    });
                }
            }
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = a_MF->const_arrays();
            auto const& mask = a_mask.const_arrays();
            mx = ParReduce(TypeList<ReduceOpMax>{}, TypeList<Real>{}, *a_MF, IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                if (!mask[box_no](i,j,k)) {
                    return AMREX_REAL_LOWEST;
                } else {
                    return ma[box_no](i,j,k,comp);
                }    
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:mx)
#endif
            for (MFIter mfi(*a_MF,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& a = a_MF->const_array(mfi);
                auto const& mask = a_mask.const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    if (mask(i,j,k)) {
                       mx = std::max(mx, a(i,j,k,comp));
                    }
                });
            }
        }
    }

    return mx;
}

// MultiFab min, exlucing EB-covered/fine-covered cells, local
Real
PeleLM::MFmin(const MultiFab *a_MF,
              const iMultiFab &a_mask,
              int comp)
{
    Real mn = std::numeric_limits<Real>::max();

#ifdef AMREX_USE_EB
    if ( a_MF->hasEBFabFactory() )
    {    
        const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_MF->Factory());
        auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& flagsma = flags.const_arrays();
            auto const& ma = a_MF->const_arrays();
            auto const& mask = a_mask.const_arrays();
            mn = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, *a_MF, IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {    
                if (flagsma[box_no](i,j,k).isCovered() ||
                    !mask[box_no](i,j,k)) {
                    return AMREX_REAL_MAX;
                } else {
                    return ma[box_no](i,j,k,comp);
                }    
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:mn)
#endif
            for (MFIter mfi(*a_MF,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                if (flags[mfi].getType(bx) != FabType::covered) {
                    auto const& flag = flags.const_array(mfi);
                    auto const& a = a_MF->const_array(mfi);
                    auto const& mask = a_mask.const_array(mfi);
                    AMREX_LOOP_3D(bx, i, j, k,
                    {
                        if (!flag(i,j,k).isCovered() && mask(i,j,k)) {
                            mn = std::min(mn, a(i,j,k,comp));
                        }
                    });
                }
            }
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_GPU
        if (Gpu::inLaunchRegion()) {
            auto const& ma = a_MF->const_arrays();
            auto const& mask = a_mask.const_arrays();
            mn = ParReduce(TypeList<ReduceOpMin>{}, TypeList<Real>{}, *a_MF, IntVect(0),
            [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<Real>
            {
                if (!mask[box_no](i,j,k)) {
                    return AMREX_REAL_MAX;
                } else {
                    return ma[box_no](i,j,k,comp);
                }    
            });
        } else
#endif
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:mn)
#endif
            for (MFIter mfi(*a_MF,true); mfi.isValid(); ++mfi) {
                Box const& bx = mfi.tilebox();
                auto const& a = a_MF->const_array(mfi);
                auto const& mask = a_mask.const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    if (mask(i,j,k)) {
                       mn = std::min(mn, a(i,j,k,comp));
                    }
                });
            }
        }
    }

    return mn;
}

// MultiLevel max, exlucing EB-covered/fine-covered cells
Vector<Real>
PeleLM::MLmax(const Vector<const MultiFab*> &a_MF,
      int scomp, int ncomp)
{
    AMREX_ASSERT(a_MF[0]->nComp() >= scomp+ncomp);

    Vector<Real> nmax(ncomp,AMREX_REAL_LOWEST);

    for (int lev = 0; lev < a_MF.size(); ++lev) {
       if (lev != finest_level) {
          for (int n = 0; n < ncomp; n++) {
             nmax[n] = std::max(nmax[n], MFmax(a_MF[lev],*m_coveredMask[lev],scomp+n));
          }
       } else {
          for (int n = 0; n < ncomp; n++) {
             nmax[n] = std::max(nmax[n], a_MF[lev]->max(scomp+n,0,true));
          }
       }
    }

    ParallelDescriptor::ReduceRealMax(nmax.data(),ncomp);
    return nmax;
}

// MultiLevel min, exlucing EB-covered/fine-covered cells
Vector<Real>
PeleLM::MLmin(const Vector<const MultiFab*> &a_MF,
      int scomp, int ncomp)
{
    AMREX_ASSERT(a_MF[0]->nComp() >= scomp+ncomp);

    Vector<Real> nmin(ncomp,AMREX_REAL_MAX);

    for (int lev = 0; lev < a_MF.size(); ++lev) {
       if (lev != finest_level) {
          for (int n = 0; n < ncomp; n++) {
             nmin[n] = std::min(nmin[n], MFmin(a_MF[lev],*m_coveredMask[lev],scomp+n));
          }
       } else {
          for (int n = 0; n < ncomp; n++) {
             nmin[n] = std::min(nmin[n], a_MF[lev]->min(scomp+n,0,true));
          }
       }
    }

    ParallelDescriptor::ReduceRealMin(nmin.data(),ncomp);
    return nmin;
}

void
PeleLM::checkMemory(const std::string &a_message)
{
    if (!m_checkMem) return;

    const int IOProc = ParallelDescriptor::IOProcessorNumber();
#ifdef AMREX_USE_GPU
    Long free_mem_avail = Gpu::Device::freeMemAvailable() / (1024*1024);
    ParallelDescriptor::ReduceLongMin(free_mem_avail, IOProc);
    Print() << "     [" << a_message << "] GPU mem. avail. (MB) " << free_mem_avail << "\n";
#else
    // MultiFab memory usage
    Long max_fab_megabytes = amrex::TotalBytesAllocatedInFabsHWM() / (1024*1024);
    ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);
    Print() << "     [" << a_message << "] MFs mem. allocated (MB) " << max_fab_megabytes << "\n";
#endif
}

#ifdef AMREX_USE_EB
// Extend the cell-centered based signed distance function
void PeleLM::extendSignedDistance(MultiFab *a_signDist,
                                  Real a_extendFactor)
{
      // This is a not-so-pretty piece of code that'll take AMReX cell-averaged
      // signed distance and propagates it manually up to the point where we need to have it
      // for derefining.
      const auto geomdata = geom[0].data();
      Real maxSignedDist = a_signDist->max(0);
      const auto& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(a_signDist->Factory());
      const auto& flags = ebfactory.getMultiEBCellFlagFab();
      int nGrowFac = flags.nGrow()+1;
       
      // First set the region far away at the max value we need
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*a_signDist,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.growntilebox();
         auto const& sd_cc = a_signDist->array(mfi);
         ParallelFor(bx, [=]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            if ( sd_cc(i,j,k) >= maxSignedDist - 1e-12 ) {
               const Real* dx = geomdata.CellSize(); 
               sd_cc(i,j,k) = nGrowFac*dx[0]*a_extendFactor;
            }
         });
      }

      // Iteratively compute the distance function in boxes, propagating accross boxes using ghost cells
      // If needed, increase the number of loop to extend the reach of the distance function
      int nMaxLoop = 4;
      for (int dloop = 1; dloop <= nMaxLoop; dloop++ ) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(*a_signDist,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            const Box& bx = mfi.tilebox();
            const Box& gbx = grow(bx,1);
            if ( flags[mfi].getType(gbx) == FabType::covered ) {
                continue;
            }
            auto const& sd_cc = a_signDist->array(mfi);
            ParallelFor(bx, [=]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               const auto glo = amrex::lbound(gbx);
               const auto ghi = amrex::ubound(gbx);
               const Real* dx = geomdata.CellSize();
               Real extendedDist = dx[0] * a_extendFactor;
               if ( sd_cc(i,j,k) >= maxSignedDist - 1e-12 ) {
                  Real closestEBDist = 1e12;
                  for (int kk = glo.z; kk <= ghi.z; ++kk) {
                     for (int jj = glo.y; jj <= ghi.y; ++jj) {
                        for (int ii = glo.x; ii <= ghi.x; ++ii) {
                           if ( (i != ii) || (j != jj) || (k != kk) ) {
                              if ( sd_cc(ii,jj,kk) > 0.0) {
                                 Real distToCell = std::sqrt( AMREX_D_TERM(  ((i-ii)*dx[0]*(i-ii)*dx[0]),
                                                                           + ((j-jj)*dx[1]*(j-jj)*dx[1]),
                                                                           + ((k-kk)*dx[2]*(k-kk)*dx[2])));
                                 Real distToEB = distToCell + sd_cc(ii,jj,kk);
                                 if ( distToEB < closestEBDist ) closestEBDist = distToEB;
                              }
                           }
                        }
                     }
                  }
                  if ( closestEBDist < 1e10 ) {
                     sd_cc(i,j,k) = closestEBDist;
                  } else {
                     sd_cc(i,j,k) = extendedDist;
                  }
               }
            });
         }
         a_signDist->FillBoundary(geom[0].periodicity());  
      }
}
#endif
