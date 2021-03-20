#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;


void PeleLM::fluxDivergence(Vector<MultiFab> &a_divergence,
                            int div_comp,
                            const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                            int flux_comp,
                            int ncomp,
                            int intensiveFluxes,
                            Real scale) {

   BL_PROFILE_VAR("PeleLM::fluxDivergence()", fluxDivergence);
   if (intensiveFluxes) {        // Fluxes are intensive -> need area scaling in div
      for (int lev = 0; lev <= finest_level; ++lev) {
         intFluxDivergenceLevel(lev,a_divergence[lev], div_comp, a_fluxes[lev], flux_comp,
                                ncomp, scale);
      }
   } else {                      // Fluxes are extensive
      for (int lev = 0; lev <= finest_level; ++lev) {
         extFluxDivergenceLevel(lev,a_divergence[lev], div_comp, a_fluxes[lev], flux_comp,
                                ncomp, scale);
      }
   }
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
      auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
#endif


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
               extFluxDivergence_K( i, j, k, ncomp,
                                    AMREX_D_DECL(fluxX, fluxY, fluxZ),
                                    vol, scale, divergence);
            } else {
               Real vfracinv = 1.0/vfrac(i,j,k);
               extFluxDivergence_K( i, j, k, ncomp,
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
            extFluxDivergence_K( i, j, k, ncomp,
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
      auto const& ebfactory = dynamic_cast<EBFArrayBoxFactory const&>(Factory());
      Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
      areafrac  = ebfactory.getAreaFrac();
#endif


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
               for (int n = 0; n < nComp; n++) {
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
               for (int n = 0; n < nComp; n++) {
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
