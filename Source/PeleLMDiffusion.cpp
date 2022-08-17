#include <PeleLM.H>
#include <PeleLMUtils.H>
#include <PeleLM_K.H>
#include <DiffusionOp.H>
#include <AMReX_PlotFileUtil.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB_utils.H>
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

DiffusionOp*
PeleLM::getDiffusionOp()
{
   if (!m_diffusion_op) m_diffusion_op.reset(new DiffusionOp(this));
   return m_diffusion_op.get();
}

DiffusionOp*
PeleLM::getMCDiffusionOp(int ncomp)
{
   if (!m_mcdiffusion_op || m_mcdiffusion_op->m_ncomp != ncomp) m_mcdiffusion_op.reset(new DiffusionOp(this,ncomp));
   return m_mcdiffusion_op.get();
}

DiffusionTensorOp*
PeleLM::getDiffusionTensorOp ()
{
    if (!m_diffusionTensor_op) m_diffusionTensor_op.reset(new DiffusionTensorOp(this));
    return m_diffusionTensor_op.get();
}

void PeleLM::computeDifferentialDiffusionTerms(const TimeStamp &a_time,
                                               std::unique_ptr<AdvanceDiffData> &diffData,
                                               int is_init)
{
   BL_PROFILE("PeleLM::computeDifferentialDiffusionTerms()");

   //----------------------------------------------------------------
   // Checks
   AMREX_ASSERT((a_time == AmrOldTime) || (a_time == AmrNewTime) );
   if (a_time == AmrOldTime) {
      AMREX_ASSERT(diffData->Dn.size() == finest_level+1);
      AMREX_ASSERT(diffData->Dn[0].nComp() >= NUM_SPECIES+2);
   } else {
      AMREX_ASSERT(diffData->Dnp1.size() == finest_level+1);
      AMREX_ASSERT(diffData->Dnp1[0].nComp() >= NUM_SPECIES+2);
   }

   //----------------------------------------------------------------
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

   //----------------------------------------------------------------
   // Compute differential diffusion fluxes including correction velocity and wbar term
   // During initialization, don't bother getting the wbar fluxes separately
   if (is_init || !m_use_wbar) {
      computeDifferentialDiffusionFluxes(a_time, GetVecOfArrOfPtrs(fluxes), {});
   } else {
      computeDifferentialDiffusionFluxes(a_time, GetVecOfArrOfPtrs(fluxes), GetVecOfArrOfPtrs(diffData->wbar_fluxes));
   }

   // If doing species balances, compute face domain integrals
   // using level 0 since we've averaged down the fluxes already
   // Factor for SDC is 0.5 is for Dn and -0.5 for Dnp1
   if ((m_sdcIter == 0 || m_sdcIter == m_nSDCmax)  && m_do_speciesBalance) {
       Real sdc_weight = (a_time == AmrOldTime) ? 0.5 : -0.5;
       addRhoYFluxes(GetArrOfConstPtrs(fluxes[0]),geom[0], sdc_weight);
   }

   //----------------------------------------------------------------
   // TODO simplify the following ...
   // Compute divergence/fill a_viscTerm
   // [0:NUM_SPECIES-1] Species           : \nabla \cdot \Flux_k
   // [NUM_SPECIES]     Temperature       : \nabla \cdot (-\lambda \nabla T)
   // [NUM_SPECIES+1]   Differential diff : \nabla \cdot \sum_k ( h_k * \Flux_k )
   if (a_time == AmrOldTime) {
#ifdef AMREX_USE_EB
      auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);
      auto bcRecSpec_d = convertToDeviceVector(bcRecSpec);
      fluxDivergenceRD(GetVecOfConstPtrs(getSpeciesVect(AmrOldTime)), 0,
                       GetVecOfPtrs(diffData->Dn), 0,
                       GetVecOfArrOfPtrs(fluxes), 0,
                       NUM_SPECIES, 1, bcRecSpec_d.dataPtr(), -1.0, m_dt);
      auto bcRecTemp = fetchBCRecArray(TEMP,1);
      auto bcRecTemp_d = convertToDeviceVector(bcRecTemp);
      fluxDivergenceRD(GetVecOfConstPtrs(getTempVect(AmrOldTime)), 0,
                       GetVecOfPtrs(diffData->Dn), NUM_SPECIES,
                       GetVecOfArrOfPtrs(fluxes), NUM_SPECIES,
                       1, 1, bcRecTemp_d.dataPtr(), -1.0, m_dt);
      auto bcRecRhoH = fetchBCRecArray(RHOH,1);
      auto bcRecRhoH_d = convertToDeviceVector(bcRecRhoH);
      fluxDivergenceRD(GetVecOfConstPtrs(getRhoHVect(AmrOldTime)), 0,
                       GetVecOfPtrs(diffData->Dn), NUM_SPECIES+1,
                       GetVecOfArrOfPtrs(fluxes), NUM_SPECIES+1,
                       1, 1, bcRecRhoH_d.dataPtr(), -1.0, m_dt);
#else
      fluxDivergence(GetVecOfPtrs(diffData->Dn), 0, GetVecOfArrOfPtrs(fluxes), 0, NUM_SPECIES+2, 1, -1.0);
#endif
   } else {
#ifdef AMREX_USE_EB
      auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);
      auto bcRecSpec_d = convertToDeviceVector(bcRecSpec);
      fluxDivergenceRD(GetVecOfConstPtrs(getSpeciesVect(AmrNewTime)), 0,
                       GetVecOfPtrs(diffData->Dnp1), 0,
                       GetVecOfArrOfPtrs(fluxes), 0,
                       NUM_SPECIES, 1, bcRecSpec_d.dataPtr(), -1.0, m_dt);
      auto bcRecTemp = fetchBCRecArray(TEMP,1);
      auto bcRecTemp_d = convertToDeviceVector(bcRecTemp);
      fluxDivergenceRD(GetVecOfConstPtrs(getTempVect(AmrNewTime)), 0,
                       GetVecOfPtrs(diffData->Dnp1), NUM_SPECIES,
                       GetVecOfArrOfPtrs(fluxes), NUM_SPECIES,
                       1, 1, bcRecTemp_d.dataPtr(), -1.0, m_dt);
      auto bcRecRhoH = fetchBCRecArray(RHOH,1);
      auto bcRecRhoH_d = convertToDeviceVector(bcRecRhoH);
      fluxDivergenceRD(GetVecOfConstPtrs(getRhoHVect(AmrNewTime)), 0,
                       GetVecOfPtrs(diffData->Dnp1), NUM_SPECIES+1,
                       GetVecOfArrOfPtrs(fluxes), NUM_SPECIES+1,
                       1, 1, bcRecRhoH_d.dataPtr(), -1.0, m_dt);
#else
      fluxDivergence(GetVecOfPtrs(diffData->Dnp1), 0, GetVecOfArrOfPtrs(fluxes), 0, NUM_SPECIES+2, 1, -1.0);
#endif
   }

   // Get the wbar term if appropriate
   if (!is_init && m_use_wbar) {
#ifdef AMREX_USE_EB
      auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);
      auto bcRecSpec_d = convertToDeviceVector(bcRecSpec);
      fluxDivergenceRD(GetVecOfConstPtrs(getSpeciesVect(a_time)), 0,
                       GetVecOfPtrs(diffData->Dwbar), 0,
                       GetVecOfArrOfPtrs(diffData->wbar_fluxes), 0,
                       NUM_SPECIES, 1, bcRecSpec_d.dataPtr(), -1.0, m_dt);
#else
      fluxDivergence(GetVecOfPtrs(diffData->Dwbar), 0, GetVecOfArrOfPtrs(diffData->wbar_fluxes), 0, NUM_SPECIES, 1, -1.0);
#endif
   }

#ifdef AMREX_USE_EB
    // Set EB-covered diffusion terms here to avoid having to mask operations using this data
    // later on
    for (int lev = 0; lev <= finest_level; ++lev) {
        if (a_time == AmrOldTime) {
            EB_set_covered(diffData->Dn[lev],0.0);
        } else {
            EB_set_covered(diffData->Dnp1[lev],0.0);
        }
        if (!is_init && m_use_wbar) {
            EB_set_covered(diffData->Dwbar[lev],0.0);
        }
    }
#endif
}

void PeleLM::computeDifferentialDiffusionFluxes(const TimeStamp &a_time,
                                                const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                                                const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_wbarfluxes)
{
   BL_PROFILE("PeleLM::computeDifferentialDiffusionFluxes()");

   //----------------------------------------------------------------
   // Species fluxes
   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

#ifdef PELE_USE_EFIELD
   // Get the species diffusion fluxes from the DiffusionOp
   // Don't average down just yet
   int do_avgDown = 0;
   getMCDiffusionOp(NUM_SPECIES-NUM_IONS)->computeDiffFluxes(a_fluxes, 0,
                                                             GetVecOfConstPtrs(getSpeciesVect(a_time)), 0,
                                                             GetVecOfConstPtrs(getDensityVect(a_time)),
                                                             GetVecOfConstPtrs(getDiffusivityVect(a_time)), 0, bcRecSpec,
                                                             NUM_SPECIES-NUM_IONS, -1.0, do_avgDown);
   // Ions one by one
   for ( int n = 0; n < NUM_IONS; n++) {
      auto bcRecIons = fetchBCRecArray(FIRSTSPEC+NUM_SPECIES-NUM_IONS+n,1);
      getDiffusionOp()->computeDiffFluxes(a_fluxes, NUM_SPECIES-NUM_IONS+n,
                                          GetVecOfConstPtrs(getSpeciesVect(a_time)), NUM_SPECIES-NUM_IONS+n,
                                          GetVecOfConstPtrs(getDensityVect(a_time)),
                                          GetVecOfConstPtrs(getDiffusivityVect(a_time)), 0, bcRecIons,
                                          1, -1.0, do_avgDown);
   }
#else
   // Get the species diffusion fluxes from the DiffusionOp
   // Don't average down just yet
   int do_avgDown = 0;
   getMCDiffusionOp(NUM_SPECIES)->computeDiffFluxes(a_fluxes, 0,
                                                    GetVecOfConstPtrs(getSpeciesVect(a_time)), 0,
                                                    GetVecOfConstPtrs(getDensityVect(a_time)),
                                                    GetVecOfConstPtrs(getDiffusivityVect(a_time)), 0, bcRecSpec,
                                                    NUM_SPECIES, -1.0, do_avgDown);
#endif

   // Add the wbar term
   if (m_use_wbar) {
      int need_wbar_fluxes = (a_wbarfluxes.empty()) ? 0 : 1;
      if ( !need_wbar_fluxes ) {
         addWbarTerm(a_fluxes,
                     {},
                     GetVecOfConstPtrs(getSpeciesVect(a_time)),
                     GetVecOfConstPtrs(getDensityVect(a_time)),
                     GetVecOfConstPtrs(getDiffusivityVect(a_time)));
      } else {
         addWbarTerm(a_fluxes,
                     a_wbarfluxes,
                     GetVecOfConstPtrs(getSpeciesVect(a_time)),
                     GetVecOfConstPtrs(getDensityVect(a_time)),
                     GetVecOfConstPtrs(getDiffusivityVect(a_time)));
      }
   }

   // Adjust species diffusion fluxes to ensure their sum is zero
   adjustSpeciesFluxes(a_fluxes,
                       GetVecOfConstPtrs(getSpeciesVect(a_time)));
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Enthalpy fluxes
   // Get the temperature BCRec
   auto bcRecTemp = fetchBCRecArray(TEMP,1);

   // Fourier: - \lambda \nabla T
   do_avgDown = 0;
   getDiffusionOp()->computeDiffFluxes(a_fluxes, NUM_SPECIES,
                                       GetVecOfConstPtrs(getTempVect(a_time)), 0,
                                       {},
                                       GetVecOfConstPtrs(getDiffusivityVect(a_time)), NUM_SPECIES, bcRecTemp,
                                       1, -1.0, do_avgDown);

   // Differential diffusion term: \sum_k ( h_k * \Flux_k )
   computeSpeciesEnthalpyFlux(a_fluxes,
                              GetVecOfConstPtrs(getTempVect(a_time)));
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Get fluxes consistent accross levels by averaging down all components
   getDiffusionOp()->avgDownFluxes(a_fluxes,0,NUM_SPECIES+2);
   //----------------------------------------------------------------

#ifdef AMREX_USE_EB
   //----------------------------------------------------------------
   // Set covered faces to large dummy values to catch any usage
   for (int lev = 0; lev <= finest_level; ++lev) {
      EB_set_covered_faces({AMREX_D_DECL(a_fluxes[lev][0],a_fluxes[lev][1],a_fluxes[lev][2])},1.234e40);
   }
#endif
}

void PeleLM::addWbarTerm(const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_spfluxes,
                         const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_spwbarfluxes,
                         Vector<MultiFab const*> const &a_spec,
                         Vector<MultiFab const*> const &a_rho,
                         Vector<MultiFab const*> const &a_beta)
{
   //------------------------------------------------------------------------
   // if a container for wbar fluxes is provided, fill it
   int need_wbar_fluxes = (a_spwbarfluxes.empty()) ? 0 : 1;

   //------------------------------------------------------------------------
   // Compute Wbar on all the levels
   int nGrow = 1;                // Need one ghost cell to compute gradWbar
   Vector<MultiFab> Wbar(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {

      Wbar[lev].define(grids[lev],dmap[lev],1,nGrow,MFInfo(),Factory(lev));

#ifdef AMREX_USE_OMP
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
   // Compute Wbar gradients and do average down to get gradients consistent accross levels
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
   getDiffusionOp()->computeGradient(GetVecOfArrOfPtrs(gradWbar),
                                     {},        // Don't need the laplacian out
                                     GetVecOfConstPtrs(Wbar),
                                     bcRecSpec[0], do_avgDown);

   //------------------------------------------------------------------------
   // add Wbar term to species fluxes
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get edge diffusivity
      int doZeroVisc = 1;
      Array<MultiFab,AMREX_SPACEDIM> beta_ec = getDiffusivity(lev, 0, NUM_SPECIES, doZeroVisc, bcRecSpec, *a_beta[lev]);

      const Box& domain = geom[lev].Domain();
      bool use_harmonic_avg = m_harm_avg_cen2edge ? true : false;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
         for (MFIter mfi(*a_beta[lev],TilingIfNotGPU()); mfi.isValid();++mfi)
         {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {

               // Get edge centered rhoYs
               const Box ebx = mfi.nodaltilebox(idim);
               FArrayBox rhoY_ed(ebx, NUM_SPECIES, The_Async_Arena());

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

               auto const& rhoY          = rhoY_ed.const_array(0);
               auto const& gradWbar_ar   = gradWbar[lev][idim].const_array(mfi);
               auto const& beta_ar       = beta_ec[idim].const_array(mfi);
               auto const& spFlux_ar     = a_spfluxes[lev][idim]->array(mfi);
               auto const& spwbarFlux_ar = ( need_wbar_fluxes ) ? a_spwbarfluxes[lev][idim]->array(mfi)
                                                                : a_spfluxes[lev][idim]->array(mfi);     // Dummy unused Array4

               // Wbar flux is : - \rho Y_m / \overline{W} * D_m * \nabla \overline{W}
               // with beta_m = \rho * D_m below
               amrex::ParallelFor(ebx, [need_wbar_fluxes, gradWbar_ar, beta_ar, rhoY, spFlux_ar, spwbarFlux_ar]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  auto eos = pele::physics::PhysicsType::eos();
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
                  eos.Y2WBAR(y, WBAR);
                  WBAR *= 0.001;
                  for (int n = 0; n < NUM_SPECIES; n++) {
                     spFlux_ar(i,j,k,n) -= y[n] / WBAR * beta_ar(i,j,k,n) * gradWbar_ar(i,j,k);
                  }
                  if ( need_wbar_fluxes ) {
                     for (int n = 0; n < NUM_SPECIES; n++) {
                        spwbarFlux_ar(i,j,k,n) = -y[n] / WBAR * beta_ar(i,j,k,n) * gradWbar_ar(i,j,k);
                     }
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

   BL_PROFILE("PeleLM::adjustSpeciesFluxes()");

   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

   for (int lev = 0; lev <= finest_level; ++lev) {

      const Box& domain = geom[lev].Domain();

#ifdef AMREX_USE_EB
      //------------------------------------------------------------------------
      // Get the edge species state needed for EB
      int nGrow = 1;
      const auto& ba = a_spec[lev]->boxArray();
      const auto& dm = a_spec[lev]->DistributionMap();
      const auto& ebfact = EBFactory(lev);
      Array<MultiFab,AMREX_SPACEDIM> edgstate{AMREX_D_DECL(MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(0)),
                                                                    dm, NUM_SPECIES, nGrow, MFInfo(), ebfact),
                                                           MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(1)),
                                                                    dm, NUM_SPECIES, nGrow, MFInfo(), ebfact),
                                                           MultiFab(amrex::convert(ba,IntVect::TheDimensionVector(2)),
                                                                    dm, NUM_SPECIES, nGrow, MFInfo(), ebfact))};
      EB_interp_CellCentroid_to_FaceCentroid(*a_spec[lev], GetArrOfPtrs(edgstate), 0, 0, NUM_SPECIES, geom[lev], bcRecSpec);
      auto const &areafrac = ebfact.getAreaFrac();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(*a_spec[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            const Box& ebx = mfi.nodaltilebox(idim);
            const Box& edomain   = amrex::surroundingNodes(domain,idim);
            auto const& rhoY     = a_spec[lev]->const_array(mfi);
            auto const& flux_dir = a_spfluxes[lev][idim]->array(mfi);

            const auto bc_lo = bcRecSpec[0].lo(idim);
            const auto bc_hi = bcRecSpec[0].hi(idim);

#ifdef AMREX_USE_EB
            auto const& flagfab  = ebfact.getMultiEBCellFlagFab()[mfi];

            if (flagfab.getType(amrex::grow(ebx,0)) != FabType::covered )
            {
               // No cut cells in tile + nghost-cell witdh halo -> use non-eb routine
               if (flagfab.getType(amrex::grow(ebx,nGrow)) == FabType::regular )
               {
                  amrex::ParallelFor(ebx, [idim, rhoY, flux_dir, edomain, bc_lo, bc_hi]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     int idx[3] = {i,j,k};
                     bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] <= edomain.smallEnd(idim) ) );
                     bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] >= edomain.bigEnd(idim) ) );
                     repair_flux( i, j, k, idim, on_lo, on_hi, rhoY, flux_dir );
                  });
               }
               else
               {
                  auto const& rhoYed_ar   = edgstate[idim].const_array(mfi);
                  auto const &areafrac_ar = areafrac[idim]->const_array(mfi);
                  amrex::ParallelFor(ebx, [idim, rhoY, flux_dir, rhoYed_ar, areafrac_ar, edomain, bc_lo, bc_hi]
                  AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                  {
                     int idx[3] = {i,j,k};
                     bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] <= edomain.smallEnd(idim) ) );
                     bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) &&
                                    ( idx[idim] >= edomain.bigEnd(idim) ) );
                     repair_flux_eb( i, j, k, idim, on_lo, on_hi, rhoY, rhoYed_ar, areafrac_ar, flux_dir );
                  });
               }
            }
#else
            amrex::ParallelFor(ebx, [idim, rhoY, flux_dir, edomain, bc_lo, bc_hi]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               int idx[3] = {i,j,k};
               bool on_lo = ( ( bc_lo == amrex::BCType::ext_dir ) &&
                              ( idx[idim] <= edomain.smallEnd(idim) ) );
               bool on_hi = ( ( bc_hi == amrex::BCType::ext_dir ) &&
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

   BL_PROFILE("PeleLM::computeSpeciesEnthalpyFlux()");

   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

   for (int lev = 0; lev <= finest_level; ++lev) {

#ifdef AMREX_USE_EB
      auto const& ebfact = EBFactory(lev);
#endif
      //------------------------------------------------------------------------
      // Compute the cell-centered species enthalpies
      int nGrow = 1;
      MultiFab Enth(grids[lev],dmap[lev],NUM_SPECIES,nGrow,MFInfo(),Factory(lev));

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Enth,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& gbx  = mfi.growntilebox();
         auto const& Temp_arr  = a_temp[lev]->const_array(mfi);
         auto const& Hi_arr    = Enth.array(mfi);

#ifdef AMREX_USE_EB
         auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
         auto const& flag    = flagfab.const_array();
         if (flagfab.getType(gbx) == FabType::covered) {              // Covered boxes
            amrex::ParallelFor(gbx, [Hi_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               Hi_arr(i,j,k) = 0.0;
            });
         } else if (flagfab.getType(gbx) != FabType::regular ) {     // EB containing boxes
            amrex::ParallelFor(gbx, [Temp_arr, Hi_arr, flag]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               if ( flag(i,j,k).isCovered() ) {
                  Hi_arr(i,j,k) = 0.0;
               } else {
                  getHGivenT( i, j, k, Temp_arr, Hi_arr );
               }
            });
         } else
#endif
         {
            amrex::ParallelFor(gbx, [Temp_arr, Hi_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               getHGivenT( i, j, k, Temp_arr, Hi_arr );
            });
         }
      }

      //------------------------------------------------------------------------
      // Get the face-centered species enthalpies
      int doZeroVisc = 0;
      Array<MultiFab,AMREX_SPACEDIM> Enth_ec = getDiffusivity(lev, 0, NUM_SPECIES, doZeroVisc, bcRecSpec, Enth);

      //------------------------------------------------------------------------
      // Compute \sum_k { \Flux_k * h_k }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(Enth,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
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
      }
   }
}

void PeleLM::differentialDiffusionUpdate(std::unique_ptr<AdvanceAdvData> &advData,
                                         std::unique_ptr<AdvanceDiffData> &diffData)
{
   BL_PROFILE("PeleLM::differentialDiffusionUpdate()");

   //------------------------------------------------------------------------
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

   //------------------------------------------------------------------------
   // Convert species forcing into actual solve RHS by *dt and adding rhoY^{n}
   // Could have done it at the same time the forcing is built, but this is clearer
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get t^{n} data pointer
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(advData->Forcing[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rhoY_o  = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& fY      = advData->Forcing[lev].array(mfi,0);
         amrex::ParallelFor(bx, [rhoY_o, fY, dt=m_dt]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            for (int n = 0; n < NUM_SPECIES; n++) {
               fY(i,j,k,n) *= dt;
               fY(i,j,k,n) += rhoY_o(i,j,k,n);
            }
         });
      }
   }

   //------------------------------------------------------------------------
   // Species diffusion solve
   // Get the species BCRec
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);

#ifdef PELE_USE_EFIELD
   // Solve for \widetilda{rhoY^{np1,kp1}}
   // -> return the uncorrected fluxes^{np1,kp1}
   // -> and the partially updated species (not including wbar or flux correction)
   getMCDiffusionOp(NUM_SPECIES-NUM_IONS)->diffuse_scalar(GetVecOfPtrs(getSpeciesVect(AmrNewTime)), 0,
                                                          GetVecOfConstPtrs(advData->Forcing), 0,
                                                          GetVecOfArrOfPtrs(fluxes), 0,
                                                          GetVecOfConstPtrs(getDensityVect(AmrNewTime)),        // this is the acoeff of LinOp
                                                          GetVecOfConstPtrs(getDensityVect(AmrNewTime)),        // this triggers proper scaling by density
                                                          GetVecOfConstPtrs(getDiffusivityVect(AmrNewTime)), 0, bcRecSpec,
                                                          NUM_SPECIES-NUM_IONS, 0, m_dt);
   // Ions one by one
   for ( int n = 0; n < NUM_IONS; n++) {
      auto bcRecIons = fetchBCRecArray(FIRSTSPEC+NUM_SPECIES-NUM_IONS+n,1);
      getDiffusionOp()->diffuse_scalar(GetVecOfPtrs(getSpeciesVect(AmrNewTime)), NUM_SPECIES-NUM_IONS+n,
                                       GetVecOfConstPtrs(advData->Forcing), NUM_SPECIES-NUM_IONS+n,
                                       GetVecOfArrOfPtrs(fluxes), NUM_SPECIES-NUM_IONS+n,
                                       GetVecOfConstPtrs(getDensityVect(AmrNewTime)),        // this is the acoeff of LinOp
                                       GetVecOfConstPtrs(getDensityVect(AmrNewTime)),        // this triggers proper scaling by density
                                       GetVecOfConstPtrs(getDiffusivityVect(AmrNewTime)), 0, bcRecIons,
                                       1, 0, m_dt);
   }
#else
   // Solve for \widetilda{rhoY^{np1,kp1}}
   // -> return the uncorrected fluxes^{np1,kp1}
   // -> and the partially updated species (not including wbar or flux correction)
   getMCDiffusionOp(NUM_SPECIES)->diffuse_scalar(GetVecOfPtrs(getSpeciesVect(AmrNewTime)), 0,
                                                 GetVecOfConstPtrs(advData->Forcing), 0,
                                                 GetVecOfArrOfPtrs(fluxes), 0,
                                                 GetVecOfConstPtrs(getDensityVect(AmrNewTime)),        // this is the acoeff of LinOp
                                                 GetVecOfConstPtrs(getDensityVect(AmrNewTime)),        // this triggers proper scaling by density
                                                 GetVecOfConstPtrs(getDiffusivityVect(AmrNewTime)), 0, bcRecSpec,
                                                 NUM_SPECIES, 0, m_dt);
#endif

   // Add lagged Wbar term
   // Computed in computeDifferentialDiffusionTerms at t^{n} if first SDC iteration, t^{np1,k} otherwise
   if (m_use_wbar) {
      for (int lev = 0; lev <= finest_level; ++lev) {

         auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
         {
            for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
            {
               const Box& ebx = mfi.nodaltilebox(idim);
               auto const& flux_spec = fluxes[lev][idim].array(mfi);
               auto const& flux_wbar = diffData->wbar_fluxes[lev][idim].const_array(mfi);
               amrex::ParallelFor(ebx, NUM_SPECIES, [ flux_spec, flux_wbar ]
               AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
               {
                   flux_spec(i,j,k,n) += flux_wbar(i,j,k,n);
               });
            }
         }
      }
   }

   // FillPatch the new species before computing flux correction terms
   fillPatchSpecies(AmrNewTime);

   // Adjust species diffusion fluxes to ensure their sum is zero
   adjustSpeciesFluxes(GetVecOfArrOfPtrs(fluxes),
                       GetVecOfConstPtrs(getSpeciesVect(AmrNewTime)));

   // Average down fluxes^{np1,kp1}
   getDiffusionOp()->avgDownFluxes(GetVecOfArrOfPtrs(fluxes),0,NUM_SPECIES);

   // Compute diffusion term D^{np1,kp1} (or Dhat)
   fluxDivergence(GetVecOfPtrs(diffData->Dhat), 0, GetVecOfArrOfPtrs(fluxes), 0, NUM_SPECIES, 1, -1.0);

   // Update species
   // Remove the Wbar term because we included it in both the dhat and the forcing.
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx     = mfi.tilebox();
         auto const& rhoY  = ldata_p->state.array(mfi,FIRSTSPEC);
         auto const& dhat  = diffData->Dhat[lev].const_array(mfi);
         auto const& force = advData->Forcing[lev].const_array(mfi,0);
         auto const& dwbar = (m_use_wbar) ? diffData->Dwbar[lev].const_array(mfi) :
                                            diffData->Dhat[lev].const_array(mfi);            // Dummy unused Array4
         amrex::ParallelFor(bx, NUM_SPECIES, [rhoY,dhat,force,dwbar,
                                              dt=m_dt,use_wbar=m_use_wbar]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            rhoY(i,j,k,n) = force(i,j,k,n) + dt * dhat(i,j,k,n);
            if (use_wbar) {
               rhoY(i,j,k,n) -= dt * dwbar(i,j,k,n);
            }
         });
      }
   }

   // FillPatch species again before going into the enthalpy solve
   fillPatchSpecies(AmrNewTime);

   // If doing species balances, compute face domain integrals
   // using level 0 since we've averaged down the fluxes already
   if (m_sdcIter == m_nSDCmax && m_do_speciesBalance) {
       addRhoYFluxes(GetArrOfConstPtrs(fluxes[0]),geom[0]);
   }
   //------------------------------------------------------------------------

   //------------------------------------------------------------------------
   // Enthalpy iterative diffusion solve
   // Get the temperature BCRec
   auto bcRecTemp = fetchBCRecArray(TEMP,1);

   // Fourier: - \lambda \nabla T
   int do_avgDown = 0;
   getDiffusionOp()->computeDiffFluxes(GetVecOfArrOfPtrs(fluxes), NUM_SPECIES,
                                       GetVecOfConstPtrs(getTempVect(AmrNewTime)), 0,
                                       {},
                                       GetVecOfConstPtrs(getDiffusivityVect(AmrNewTime)), NUM_SPECIES, bcRecTemp,
                                       1, -1.0, do_avgDown);

   // Differential diffusion term: \sum_k ( h_k * \Flux_k )
   computeSpeciesEnthalpyFlux(GetVecOfArrOfPtrs(fluxes),
                              GetVecOfConstPtrs(getTempVect(AmrNewTime)));

   // average_down enthalpy fluxes
   getDiffusionOp()->avgDownFluxes(GetVecOfArrOfPtrs(fluxes),NUM_SPECIES,2);

   // Compute diffusion term D^{np1,kp1} of Fourier and DifferentialDiffusion
   fluxDivergence(GetVecOfPtrs(diffData->Dhat), NUM_SPECIES, GetVecOfArrOfPtrs(fluxes), NUM_SPECIES, 2, 1, -1.0);

   //------------------------------------------------------------------------
   // delta(T) iterations
   if (m_deltaT_verbose) {
      Print() << " Iterative solve for deltaT \n";
   }

   //------------------------------------------------------------------------
   // Temporary data holders
   Vector<MultiFab> rhs(finest_level+1);                 // Linear deltaT solve RHS
   Vector<MultiFab> Tsave(finest_level+1);               // Storage of T while working on deltaT
   Vector<MultiFab> RhoCp(finest_level+1);               // Acoeff of the linear solve
   for (int lev = 0; lev <= finest_level; ++lev) {
      rhs[lev].define(grids[lev],dmap[lev], 1, 0, MFInfo(), Factory(lev));
      Tsave[lev].define(grids[lev],dmap[lev], 1, 1, MFInfo(), Factory(lev));
      RhoCp[lev].define(grids[lev],dmap[lev], 1, 0, MFInfo(), Factory(lev));
   }

   // DeltaT norm
   Real deltaT_norm = 0.0;
   for (int dTiter = 0; dTiter < m_deltaTIterMax && (dTiter==0 || deltaT_norm >= m_deltaT_norm_max); ++dTiter) {

      // Prepare the deltaT iteration linear solve:
      // -> Assemble the RHS
      // -> Compute current value of \rho * \Cp_{mix}
      // -> Save current T^{np1,kp1}
      // -> set T^{np1,kp1} to zero
      deltaTIter_prepare(GetVecOfPtrs(rhs),
                         GetVecOfPtrs(Tsave),
                         GetVecOfPtrs(RhoCp),
                         advData, diffData);

      // Diffuse deltaT
      getDiffusionOp()->diffuse_scalar(GetVecOfPtrs(getTempVect(AmrNewTime)), 0,
                                       GetVecOfConstPtrs(rhs), 0,
                                       GetVecOfArrOfPtrs(fluxes), NUM_SPECIES,
                                       GetVecOfConstPtrs(RhoCp),
                                       {},
                                       GetVecOfConstPtrs(getDiffusivityVect(AmrNewTime)), NUM_SPECIES, bcRecTemp,
                                       1, 0, m_dt);

      // Post deltaT iteration linear solve
      // -> evaluate deltaT_norm
      // -> add deltaT to T^{np1,kp1}
      // -> recompute enthalpy fluxes
      // -> recompute rhoH
      deltaTIter_update(dTiter,
                        GetVecOfArrOfPtrs(fluxes),
                        GetVecOfConstPtrs(Tsave),
                        diffData, deltaT_norm);

      // Check for convergence failure
      if ( (dTiter == m_deltaTIterMax-1) && ( deltaT_norm > m_deltaT_norm_max ) ) {
         if ( m_crashOnDeltaTFail ) {
            Abort("deltaT_iters not converged !");
         } else {
            Print() << "deltaT_iters not converged !\n";
         }
      }
   }
   //------------------------------------------------------------------------
}

void PeleLM::deltaTIter_prepare(const Vector<MultiFab*> &a_rhs,
                                const Vector<MultiFab*> &a_Tsave,
                                const Vector<MultiFab*> &a_rhoCp,
                                std::unique_ptr<AdvanceAdvData> &advData,
                                std::unique_ptr<AdvanceDiffData> &diffData)
{
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNew_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx        = mfi.tilebox();
         // RHS pieces
         auto const& rhoH_o   = ldataOld_p->state.const_array(mfi,RHOH);
         auto const& rhoH_n   = ldataNew_p->state.const_array(mfi,RHOH);
         auto const& force    = advData->Forcing[lev].const_array(mfi,NUM_SPECIES);
         auto const& fourier  = diffData->Dhat[lev].const_array(mfi,NUM_SPECIES);
         auto const& diffDiff = diffData->Dhat[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& rhs      = a_rhs[lev]->array(mfi);
         Real dtinv           = 1.0/m_dt;

         // Cpmix
         auto const& rho      = ldataNew_p->state.const_array(mfi,DENSITY);
         auto const& rhoY     = ldataNew_p->state.const_array(mfi,FIRSTSPEC);
         auto const& T        = ldataNew_p->state.const_array(mfi,TEMP);
         auto const& rhocp    = a_rhoCp[lev]->array(mfi);

         // T save
         auto const& tsave    = a_Tsave[lev]->array(mfi);
         amrex::ParallelFor(bx, [=,dt=m_dt]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            // Assemble deltaT RHS
            rhs(i,j,k) = dt * (( rhoH_o(i,j,k) - rhoH_n(i,j,k) ) * dtinv
                         + force(i,j,k) + fourier(i,j,k) + diffDiff(i,j,k));

            // Get \rho * Cp_{mix}
            getCpmixGivenRYT( i, j, k, rho, rhoY, T, rhocp );
            rhocp(i,j,k) *= rho(i,j,k);

            // Save T
            tsave(i,j,k) = T(i,j,k);
         });
      }

      // Set T^{np1} to zero
      // Include one ghost cell to ensure levelBC at zero for linear solve
      ldataNew_p->state.setVal(0.0,TEMP,1,1);
   }
}

void PeleLM::deltaTIter_update(int a_dtiter,
                               const Vector<Array<MultiFab*,AMREX_SPACEDIM> > &a_fluxes,
                               const Vector<MultiFab const*> &a_Tsave,
                               std::unique_ptr<AdvanceDiffData> &diffData,
                               Real &a_deltaT_norm)
{

   //------------------------------------------------------------------------
   // Evaluate deltaT norm and add Tsave back into the new LevelData
   a_deltaT_norm = -1.0e12;
   for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      a_deltaT_norm = std::max(a_deltaT_norm,ldata_p->state.norm0(TEMP,0,false,true));
      MultiFab::Add(ldata_p->state,*a_Tsave[lev],0,TEMP,1,0);
   }

   if (m_deltaT_verbose) {
      Print() << "   DeltaT solve norm [" << a_dtiter << "] = " << a_deltaT_norm << "\n";
   }

   // FillPatch the new temperature before going into the fluxe computation
   fillPatchTemp(AmrNewTime);

   //------------------------------------------------------------------------
   // Re-evaluate the enthalpy fluxes and divergence
   // Get the temperature BCRec
   auto bcRecTemp = fetchBCRecArray(TEMP,1);

   // Fourier: - \lambda \nabla T
   int do_avgDown = 0;
   getDiffusionOp()->computeDiffFluxes(a_fluxes, NUM_SPECIES,
                                       GetVecOfConstPtrs(getTempVect(AmrNewTime)), 0,
                                       {},
                                       GetVecOfConstPtrs(getDiffusivityVect(AmrNewTime)), NUM_SPECIES, bcRecTemp,
                                       1, -1.0, do_avgDown);

   // Differential diffusion term: \sum_k ( h_k * \Flux_k )
   computeSpeciesEnthalpyFlux(a_fluxes,
                              GetVecOfConstPtrs(getTempVect(AmrNewTime)));

   // average_down enthalpy fluxes
   getDiffusionOp()->avgDownFluxes(a_fluxes,NUM_SPECIES,2);

   // Compute diffusion term D^{np1,kp1} of Fourier and DifferentialDiffusion
   fluxDivergence(GetVecOfPtrs(diffData->Dhat), NUM_SPECIES, a_fluxes, NUM_SPECIES, 2, 1, -1.0);

   //------------------------------------------------------------------------
   // Recompute RhoH
   for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      auto const& sma = ldata_p->state.arrays();
      amrex::ParallelFor(ldata_p->state, [=]
      AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
      {
         getRHmixGivenTY( i,j,k,
                          Array4<Real const>(sma[box_no],DENSITY),
                          Array4<Real const>(sma[box_no],FIRSTSPEC),
                          Array4<Real const>(sma[box_no],TEMP),
                          Array4<Real      >(sma[box_no],RHOH));
      });
   }
   Gpu::streamSynchronize();
}

void PeleLM::getScalarDiffForce(std::unique_ptr<AdvanceAdvData> &advData,
                                std::unique_ptr<AdvanceDiffData> &diffData)
{
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get t^{n} data pointer
      auto ldataR_p = getLevelDataReactPtr(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(advData->Forcing[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& dn      = diffData->Dn[lev].const_array(mfi,0);
         auto const& ddn     = diffData->Dn[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& dnp1k   = diffData->Dnp1[lev].const_array(mfi,0);
         auto const& ddnp1k  = diffData->Dnp1[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& r       = ldataR_p->I_R.const_array(mfi);
         auto const& a       = advData->AofS[lev].const_array(mfi,FIRSTSPEC);
         auto const& extRhoY = m_extSource[lev]->const_array(mfi,FIRSTSPEC);
         auto const& extRhoH = m_extSource[lev]->const_array(mfi,RHOH);
         auto const& fY      = advData->Forcing[lev].array(mfi,0);
         auto const& fT      = advData->Forcing[lev].array(mfi,NUM_SPECIES);
         auto const& dwbar   = (m_use_wbar) ? diffData->Dwbar[lev].const_array(mfi,0)
                                            : diffData->Dn[lev].const_array(mfi,0);          // Dummy unsed Array4
         amrex::ParallelFor(bx, [dn, ddn, dnp1k, ddnp1k, dwbar, use_wbar=m_use_wbar,do_react=m_do_react,
                                 r, a, extRhoY, extRhoH, fY, fT, dp0dt=m_dp0dt, is_closed_ch=m_closed_chamber]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            buildDiffusionForcing( i, j, k, dn, ddn, dnp1k, ddnp1k, r, a, dp0dt, is_closed_ch, do_react, fY, fT );
            if (use_wbar) {
               for (int n = 0; n < NUM_SPECIES; n++) {
                  fY(i,j,k,n) += dwbar(i,j,k,n);
               }
            }
            for (int n = 0; n < NUM_SPECIES; n++) {
              fY(i,j,k,n) += extRhoY(i,j,k,n);
            }
            fT(i,j,k) += extRhoH(i,j,k);
         });
      }
   }

   // Fill forcing ghost cells
   if ( advData->Forcing[0].nGrow() > 0 ) {
      fillpatch_forces(m_cur_time, GetVecOfPtrs(advData->Forcing), advData->Forcing[0].nGrow());
   }
}

void PeleLM::computeDivTau(const TimeStamp &a_time,
                           const Vector<MultiFab*> &a_divtau,
                           int use_density,
                           Real scale)
{
   BL_PROFILE("PeleLM::computeDivTau()");
   // Get the density component BCRec to get viscosity on faces
   auto bcRec = fetchBCRecArray(DENSITY,1);

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

void PeleLM::diffuseVelocity()
{
   BL_PROFILE("PeleLM::diffuseVelocity()");
   // Get the density component BCRec to get viscosity on faces
   auto bcRec = fetchBCRecArray(DENSITY,1);

   // CrankNicholson 0.5 coeff
   Real dt_lcl = 0.5 * m_dt;
   if (m_incompressible) {
      getDiffusionTensorOp()->diffuse_velocity(GetVecOfPtrs(getVelocityVect(AmrNewTime)),
                                               {},
                                               GetVecOfConstPtrs(getViscosityVect(AmrNewTime)),
                                               bcRec[0], dt_lcl);
   } else {
      getDiffusionTensorOp()->diffuse_velocity(GetVecOfPtrs(getVelocityVect(AmrNewTime)),
                                               GetVecOfConstPtrs(getDensityVect(AmrHalfTime)),
                                               GetVecOfConstPtrs(getViscosityVect(AmrNewTime)),
                                               bcRec[0], dt_lcl);
   }
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
