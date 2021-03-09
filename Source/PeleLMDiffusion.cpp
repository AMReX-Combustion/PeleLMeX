#include <PeleLM.H>
#include <PeleLM_K.H>
#include <DiffusionOp.H>

using namespace amrex;

DiffusionOp*
PeleLM::getDiffusionOp()
{
   if (!m_diffusion_op) m_diffusion_op.reset(new DiffusionOp(this));
   return m_diffusion_op.get();
}

DiffusionTensorOp*
PeleLM::getDiffusionTensorOp ()
{
    if (!m_diffusionTensor_op) m_diffusionTensor_op.reset(new DiffusionTensorOp(this));
    return m_diffusionTensor_op.get();
}

void PeleLM::computeDifferentialDiffusionTerms(const TimeStamp &a_time,
                                               Vector<MultiFab> &a_viscTerm)
{
   BL_PROFILE_VAR("PeleLM::computeDifferentialDiffusionTerms()", computeDifferentialDiffusionTerms);

   AMREX_ASSERT(a_viscTerm.size() == finest_level+1);
   AMREX_ASSERT(a_viscTerm[0].nComp() >= NUM_SPECIES+2);

   // Setup fluxes
   // [0:NUM_SPECIES-1] Species     : \Flux_k
   // [NUM_SPECIES]     Temperature : - \lambda \nabla T
   // [NUM_SPECIES+1]   DiffDiff    : \sum_k ( h_k * \Flux_k )
   int nGrow = 0;                   // No need for ghost face on fluxes
   Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      const auto& ba = grids[lev];
      const auto& factory = Factory(lev);
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         fluxes[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                  dmap[lev], NUM_SPECIES+2, nGrow, MFInfo(), factory);
      }
   }

   // Compute differential diffusion fluxes including correction velocity and wbar term
   computeDifferentialDiffusionFluxes(a_time, GetVecOfArrOfPtrs(fluxes));

   // Compute divergence/fill a_viscTerm
   // [0:NUM_SPECIES-1] Species           : \nabla \cdot \Flux_k
   // [NUM_SPECIES]     Temperature       : \nabla \cdot (-\lambda \nabla T)
   // [NUM_SPECIES+1]   Differential diff : \nabla \cdot \sum_k ( h_k * \Flux_k )
   fluxDivergence(a_viscTerm, 0, GetVecOfArrOfPtrs(fluxes), 0, NUM_SPECIES+2, -1.0);

}

void PeleLM::computeDifferentialDiffusionFluxes(const TimeStamp &a_time,
                                                const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes)
{
   BL_PROFILE_VAR("PeleLM::computeDifferentialDiffusionFluxes()", computeDifferentialDiffusionFluxes);

   //----------------------------------------------------------------
   // Species fluxes
   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

   // Get the species diffusion fluxes from the DiffusionOp
   // Don't average down just yet
   int do_avgDown = 0;
   getDiffusionOp()->computeDiffFluxes(a_fluxes, 0,
                                       GetVecOfConstPtrs(getSpeciesVect(a_time)), 0,
                                       GetVecOfConstPtrs(getDensityVect(a_time)),
                                       GetVecOfConstPtrs(getDiffusivityVect(a_time)), 0, bcRecSpec,
                                       NUM_SPECIES, 1.0, do_avgDown);

   // Add the wbar term
   if (m_use_wbar) {
      addWbarTerm(a_fluxes,
                  GetVecOfConstPtrs(getSpeciesVect(a_time)),
                  GetVecOfConstPtrs(getDensityVect(a_time)),
                  GetVecOfConstPtrs(getDiffusivityVect(a_time)));
   }

   // Adjust species diffusion fluxes to ensure their sum is zero
   adjustSpeciesFluxes(a_fluxes,
                       GetVecOfConstPtrs(getSpeciesVect(a_time)));

   // Get fluxes consistent accross levels by averaging down
   getDiffusionOp()->avgDownFluxes(a_fluxes,0,NUM_SPECIES);
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Enthalpy fluxes
   // Get the temperature BCRec
   auto bcRecTemp = fetchBCRecArray(TEMP,1);

   // Fourier: - \lambda \nabla T
   // Do the average down right now
   do_avgDown = 1;
   getDiffusionOp()->computeDiffFluxes(a_fluxes, NUM_SPECIES,
                                       GetVecOfConstPtrs(getTempVect(a_time)), 0,
                                       {},
                                       GetVecOfConstPtrs(getDiffusivityVect(a_time)), NUM_SPECIES, bcRecTemp, 
                                       1, -1.0, do_avgDown);

   // Differential diffusion term: \sum_k ( h_k * \Flux_k ) 
   computeSpeciesEnthalpyFlux(a_fluxes,
                              GetVecOfConstPtrs(getTempVect(a_time)));
   //----------------------------------------------------------------
}

void PeleLM::addWbarTerm(const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_spfluxes,
                         Vector<MultiFab const*> const &a_spec,
                         Vector<MultiFab const*> const &a_rho,
                         Vector<MultiFab const*> const &a_beta)
{

   //------------------------------------------------------------------------
   // Compute Wbar on all the levels
   int nGrow = 1;                // Need one ghost cell to compute gradWbar
   Vector<MultiFab> Wbar(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {

      Wbar[lev].define(grids[lev],dmap[lev],1,nGrow,MFInfo(),Factory(lev));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Wbar[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx = mfi.growntilebox();
         auto const& rho_arr  = a_rho[lev]->const_array(mfi);
         auto const& rhoY_arr = a_spec[lev]->const_array(mfi);
         auto const& Wbar_arr = Wbar[lev].array(mfi);
         amrex::ParallelFor(gbx, [rho_arr, rhoY_arr, Wbar_arr]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getMwmixGivenRY(i, j, k, rho_arr, rhoY_arr, Wbar_arr);
         });
      }
   }

   //------------------------------------------------------------------------
   // Compute Wbar gradients and do average down to get graidnet consistent accross levels
   // Get the species BCRec
   int do_avgDown = 1;
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

   nGrow = 0;                            // No need for ghost face on fluxes
   Vector<Array<MultiFab,AMREX_SPACEDIM> > gradWbar(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      const auto& ba = grids[lev];
      const auto& factory = Factory(lev);
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         gradWbar[lev][idim].define(amrex::convert(ba,IntVect::TheDimensionVector(idim)),
                                    dmap[lev], NUM_SPECIES, nGrow, MFInfo(), factory);
         gradWbar[lev][idim].setVal(0.0);
      }
   }
   getDiffusionOp()->computeGradient(GetVecOfArrOfPtrs(gradWbar), GetVecOfConstPtrs(Wbar), bcRecSpec[0], do_avgDown);

   //------------------------------------------------------------------------
   // add Wbar term to species fluxes
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get edge diffusivity
      Array<MultiFab,AMREX_SPACEDIM> beta_ec = getDiffusivity(lev, 0, NUM_SPECIES, bcRecSpec, *a_beta[lev]);

      const Box& domain = geom[lev].Domain();
      bool use_harmonic_avg = m_harm_avg_cen2edge ? true : false;

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
         FArrayBox rhoY_ed;
         for (MFIter mfi(*a_beta[lev],TilingIfNotGPU()); mfi.isValid();++mfi)
         {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {

               // Get edge centered rhoYs
               const Box ebx = mfi.nodaltilebox(idim);
               rhoY_ed.resize(ebx,NUM_SPECIES);
               Elixir rhoY_el = rhoY_ed.elixir();

               const Box& edomain = amrex::surroundingNodes(domain,idim);
               auto const& rhoY_arr = a_spec[lev]->const_array(mfi);
               const auto& rhoYed_arr = rhoY_ed.array(0);
               const auto bc_lo = bcRecSpec[0].lo(idim);
               const auto bc_hi = bcRecSpec[0].hi(idim);
               amrex::ParallelFor(ebx, [idim, bc_lo, bc_hi, use_harmonic_avg, rhoY_arr, rhoYed_arr, edomain]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  int idx[3] = {i,j,k};
                  bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) &&
                                 ( idx[idim] <= edomain.smallEnd(idim) ) );
                  bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) &&
                                 ( idx[idim] >= edomain.bigEnd(idim) ) );
                  cen2edg_cpp( i, j, k, idim, NUM_SPECIES, use_harmonic_avg, on_lo, on_hi, rhoY_arr, rhoYed_arr);
               });

               auto const& rhoY        = rhoY_ed.const_array(0);
               auto const& gradWbar_ar = gradWbar[lev][idim].const_array(mfi);
               auto const& beta_ar     = beta_ec[idim].const_array(mfi);
               auto const& spFlux_ar   = a_spfluxes[lev][idim]->array(mfi);

               // Wbar flux is : - \rho Y_m / \overline{W} * D_m * \nabla \overline{W}
               // with beta_m = \rho * D_m below
               amrex::ParallelFor(ebx, [gradWbar_ar, beta_ar, rhoY, spFlux_ar]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  // Get Wbar from rhoYs
                  amrex::Real rho = 0.0;
                  for (int n = 0; n < NUM_SPECIES; n++) {
                     rho += rhoY(i,j,k,n);
                  }
                  amrex::Real rho_inv = 1.0 / rho;
                  amrex::Real y[NUM_SPECIES] = {0.0};
                  for (int n = 0; n < NUM_SPECIES; n++) {
                     y[n] = rhoY(i,j,k,n) * rho_inv;
                  }
                  amrex::Real WBAR = 0.0;
                  EOS::Y2WBAR(y, WBAR);
                  WBAR *= 0.001;
                  for (int n = 0; n < NUM_SPECIES; n++) {
                     spFlux_ar(i,j,k,n) -= y[n] / WBAR * beta_ar(i,j,k,n) * gradWbar_ar(i,j,k);
                  }
               });
            }
         }
      }
   }
}

void PeleLM::adjustSpeciesFluxes(const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_spfluxes,
                                 Vector<MultiFab const*> const &a_spec)
{

   BL_PROFILE_VAR("PeleLM::adjustSpeciesFluxes()", adjustSpeciesFluxes);

   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

   for (int lev = 0; lev <= finest_level; ++lev) {

      const Box& domain = geom[lev].Domain();

#ifdef AMREX_USE_EB
      // Get the edge species state needed for EB
      int nGrow = 1;
      const auto& ba = a_spec[lev]->boxArray();
      const auto& dm = a_spec[lev]->DistributionMap();
      const auto& factory = a_spec[lev]->Factory();
      Array<MultiFab,AMREX_SPACEDIM> edgstate{AMREX_D_DECL(MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                                                    dm, NUM_SPECIES, nGrow, MFInfo(), factory),
                                                           MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                                                    dm, NUM_SPECIES, nGrow, MFInfo(), factory),
                                                           MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                                                    dm, NUM_SPECIES, nGrow, MFInfo(), factory))};
      EB_interp_CellCentroid_to_FaceCentroid(*a_spec[lev], GetArrOfPtrs(edgstate), 0, 0, NUM_SPECIES, geom[lev], bcRecSpec);
#endif

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*a_spec[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Box& ebx = mfi.nodaltilebox(idim);
            const Box& edomain   = amrex::surroundingNodes(domain,idim);
            auto const& rhoY     = a_spec[lev]->const_array(mfi);
            auto const& flux_dir = a_spfluxes[lev][idim]->array(mfi);

#ifdef AMREX_USE_EB
            const EBFArrayBox&  state_fab = static_cast<EBFArrayBox const&>(a_spec[lev][mfi]);
            const EBCellFlagFab&    flags = state_fab.getEBCellFlagFab();

            if (flags.getType(amrex::grow(ebx,0)) != FabType::covered )
            {
               // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
               if (flags.getType(amrex::grow(ebx,nghost)) == FabType::regular )
               {
                  amrex::ParallelFor(ebx, [idim, rhoY, flux_dir, edomain, bcRecSpec]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     int idx[3] = {i,j,k};
                     bool on_lo = ( ( bcRecSpec[0].lo(idim) == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] <= edomain.smallEnd(idim) ) );
                     bool on_hi = ( ( bcRecSpec[0].hi(idim) == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] >= edomain.bigEnd(idim) ) );
                     repair_flux( i, j, k, idim, on_lo, on_hi, rhoY, flux_dir );
                  });
               }
               else
               {
                  auto const& rhoYed_ar   = edgstate[idim].const_array(mfi);
                  auto const& areafrac_ar = areafrac[idim]->const_array(mfi);
                  amrex::ParallelFor(ebx, [idim, rhoY, flux_dir, rhoYed_d, areafrac_d, edomain, bcRecSpec]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     int idx[3] = {i,j,k};
                     bool on_lo = ( ( bcRecSpec[0].lo(idim) == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] <= edomain.smallEnd(idim) ) );
                     bool on_hi = ( ( bcRecSpec[0].hi(idim) == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] >= edomain.bigEnd(idim) ) );
                     repair_flux_eb( i, j, k, idim, on_lo, on_hi, rhoY, rhoYed_ar, areafrac_ar, flux_dir );
                  });
               }
            }
#else
            amrex::ParallelFor(ebx, [idim, rhoY, flux_dir, edomain, bcRecSpec]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               int idx[3] = {i,j,k};
               bool on_lo = ( ( bcRecSpec[0].lo(idim) == amrex::BCType::ext_dir ) &&
                              ( idx[idim] <= edomain.smallEnd(idim) ) );
               bool on_hi = ( ( bcRecSpec[0].hi(idim) == amrex::BCType::ext_dir ) &&
                              ( idx[idim] >= edomain.bigEnd(idim) ) );
               repair_flux( i, j, k, idim, on_lo, on_hi, rhoY, flux_dir );
            });
#endif
         }
      }
   }
}

void PeleLM::computeSpeciesEnthalpyFlux(const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                                        Vector<MultiFab const*> const &a_temp)
{

   BL_PROFILE_VAR("PeleLM::computeSpeciesEnthalpyFlux()", computeSpeciesEnthalpyFlux);

   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

   for (int lev = 0; lev <= finest_level; ++lev) {

      //------------------------------------------------------------------------
      // Compute the cell-centered species enthalpies
      int nGrow = 1;
      MultiFab Enth(grids[lev],dmap[lev],NUM_SPECIES,nGrow,MFInfo(),Factory(lev));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Enth,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx  = mfi.growntilebox();
         auto const& Temp_arr  = a_temp[lev]->const_array(mfi);
         auto const& Hi_arr    = Enth.array(mfi);
         amrex::ParallelFor(gbx, [Temp_arr, Hi_arr]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            getHGivenT( i, j, k, Temp_arr, Hi_arr );
         });
      }

      //------------------------------------------------------------------------
      // Get the face-centered species enthalpies
      Array<MultiFab,AMREX_SPACEDIM> Enth_ec = getDiffusivity(lev, 0, NUM_SPECIES, bcRecSpec, Enth);

      //------------------------------------------------------------------------
      // Compute \sum_k { \Flux_k * h_k }

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Enth,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
#if AMREX_USE_EB
         // TODO
#else
         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
           const Box&  ebox        = mfi.nodaltilebox(idim);
           auto const& spflux_ar   = a_fluxes[lev][idim]->const_array(mfi, 0);
           auto const& enthflux_ar = a_fluxes[lev][idim]->array(mfi, NUM_SPECIES+1);
           auto const& enth_ar     = Enth_ec[idim].const_array(mfi);
           amrex::ParallelFor(ebox, [spflux_ar, enthflux_ar, enth_ar]
           AMREX_GPU_DEVICE (int i, int j, int k) noexcept
           {
              enthflux_ar(i,j,k) = 0.0;
              for (int n = 0; n < NUM_SPECIES; n++) {
                  enthflux_ar(i,j,k) += spflux_ar(i,j,k,n)*enth_ar(i,j,k,n);
              }
           });
         } 
#endif
      }
   }
}

void PeleLM::computeDivTau(const TimeStamp &a_time,
                           const Vector<MultiFab*> &a_divtau,
                           int use_density,
                           Real scale)
{
   // Get the first velocity component BCRec
   auto bcRec = fetchBCRecArray(VELX,1);

   if (use_density) {
      getDiffusionTensorOp()->compute_divtau(a_divtau,
                                             GetVecOfConstPtrs(getVelocityVect(a_time)),
                                             GetVecOfConstPtrs(getDensityVect(a_time)),
                                             GetVecOfConstPtrs(getViscosityVect(a_time)),
                                             bcRec[0], scale);
   } else {
      getDiffusionTensorOp()->compute_divtau(a_divtau,
                                             GetVecOfConstPtrs(getVelocityVect(a_time)),
                                             {},
                                             GetVecOfConstPtrs(getViscosityVect(a_time)),
                                             bcRec[0], scale);
   }
}

void PeleLM::diffuseVelocity(const TimeStamp &a_time) 
{

   // Get rho halftime


}

Array<LinOpBCType,AMREX_SPACEDIM>
PeleLM::getDiffusionLinOpBC(Orientation::Side a_side, const BCRec &a_bc ) {

   Array<LinOpBCType,AMREX_SPACEDIM> r;
   for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (Geom(0).isPeriodic(idim)) {
         r[idim] = LinOpBCType::Periodic;
      } else {
         auto amrexbc = (a_side==Orientation::low) ? a_bc.lo(idim) : a_bc.hi(idim);
         if (amrexbc == EXT_DIR) {
            r[idim] = LinOpBCType::Dirichlet;
         } else if ( amrexbc == FOEXTRAP ||
                     amrexbc == HOEXTRAP ||
                     amrexbc == REFLECT_EVEN ) {
            r[idim] = LinOpBCType::Neumann;
         } else if ( amrexbc == REFLECT_ODD ) {
            r[idim] = LinOpBCType::reflect_odd;
         } else {
            r[idim] = LinOpBCType::bogus;
         }
      }
   }
   return r;
}

Vector<Array<LinOpBCType,AMREX_SPACEDIM>>
PeleLM::getDiffusionTensorOpBC(Orientation::Side a_side, const Vector<BCRec> a_bc) {
   AMREX_ASSERT(a_bc.size() == AMREX_SPACEDIM);
   Vector<Array<LinOpBCType,AMREX_SPACEDIM>> r(AMREX_SPACEDIM);
   for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (Geom(0).isPeriodic(idim)) {
         AMREX_D_TERM(r[0][idim] = LinOpBCType::Periodic;,
                      r[1][idim] = LinOpBCType::Periodic;,
                      r[2][idim] = LinOpBCType::Periodic;);
      } else {
         for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            auto amrexbc = (a_side==Orientation::low) ? a_bc[dir].lo(idim) : a_bc[dir].hi(idim);
            if (amrexbc == EXT_DIR) {
               r[dir][idim] = LinOpBCType::Dirichlet;
            } else if ( amrexbc == FOEXTRAP ||
                        amrexbc == HOEXTRAP ||
                        amrexbc == REFLECT_EVEN ) {
               r[dir][idim] = LinOpBCType::Neumann;
            } else if ( amrexbc == REFLECT_ODD ) {
               r[dir][idim] = LinOpBCType::reflect_odd;
            } else {
               r[dir][idim] = LinOpBCType::bogus;
            }
         }
      }
   }
   return r;
}
