#include <PeleLM.H>
#include <PeleLM_K.H>
#include <PeleLMUtils.H>
#include <hydro_utils.H>

using namespace amrex;

void PeleLM::computeVelocityAdvTerm(std::unique_ptr<AdvanceAdvData> &advData)
{

   //----------------------------------------------------------------
   // Create temporary containers
   int nGrow_force = 1;
   Vector<MultiFab> divtau(finest_level+1);
   Vector<MultiFab> velForces(finest_level+1);
   Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
   Vector<Array<MultiFab,AMREX_SPACEDIM> > faces(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), Factory(lev));
      velForces[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force, MFInfo(), Factory(lev));
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         fluxes[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                  dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), Factory(lev));
         faces[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                 dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), Factory(lev));
      }
   }

   //----------------------------------------------------------------
   // Get viscous forces
   int use_density = 0;
   computeDivTau(AmrOldTime,GetVecOfPtrs(divtau),use_density);

   //----------------------------------------------------------------
   // Gather all the velocity forces
   // F = [ (gravity+...) - gradP + divTau ] / rho
   int add_gradP = 1;
   getVelForces(AmrOldTime,GetVecOfPtrs(divtau),GetVecOfPtrs(velForces),nGrow_force,add_gradP);

   auto bcRecVel = fetchBCRecArray(VELX,AMREX_SPACEDIM);
   auto bcRecVel_d = convertToDeviceVector(bcRecVel);
   auto AdvTypeVel = fetchAdvTypeArray(VELX,AMREX_SPACEDIM);
   auto AdvTypeVel_d = convertToDeviceVector(AdvTypeVel);
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);

      //----------------------------------------------------------------
      // Get divU
      int nGrow_divu = 4;  // Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         Real time  = getTime(lev,AmrOldTime);
         fillpatch_divu(lev,time,divu,nGrow_divu);
      }

      //----------------------------------------------------------------
#ifdef AMREX_USE_EB
      const auto& ebfact = EBFactory(lev);
#endif

      //----------------------------------------------------------------
      // Compute the velocity fluxes
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& fx = fluxes[lev][0].array(mfi);,
                      auto const& fy = fluxes[lev][1].array(mfi);,
                      auto const& fz = fluxes[lev][2].array(mfi);)
         AMREX_D_TERM(auto const& facex = faces[lev][0].array(mfi);,
                      auto const& facey = faces[lev][1].array(mfi);,
                      auto const& facez = faces[lev][2].array(mfi);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& vel_arr   = ldata_p->state.const_array(mfi,VELX);
         auto const& force_arr = velForces[lev].const_array(mfi);

         bool is_velocity = true;
         bool fluxes_are_area_weighted = false;
         bool knownEdgeState = false;
         HydroUtils::ComputeFluxesOnBoxFromState(bx, AMREX_SPACEDIM, mfi,
                                                 vel_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 AMREX_D_DECL(facex,facey,facez),
                                                 knownEdgeState,
                                                 AMREX_D_DECL(umac, vmac, wmac),
                                                 divu_arr, force_arr,
                                                 geom[lev], m_dt,
                                                 bcRecVel, bcRecVel_d.dataPtr(), AdvTypeVel_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                 ebfact,
#endif
                                                 m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                 is_velocity, fluxes_are_area_weighted,
                                                 m_advection_type);
      }
#ifdef AMREX_USE_EB
      EB_set_covered_faces(GetArrOfPtrs(fluxes[lev]),0.);
      EB_set_covered_faces(GetArrOfPtrs(faces[lev]),0.);
#endif
   }

   //----------------------------------------------------------------
   // Average down fluxes to ensure C/F consistency
   for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
      EB_average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                            GetArrOfPtrs(fluxes[lev-1]),
                            refRatio(lev-1),geom[lev-1]);
      EB_average_down_faces(GetArrOfConstPtrs(faces[lev]),
                            GetArrOfPtrs(faces[lev-1]),
                            refRatio(lev-1),geom[lev-1]);
#else
      average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                         GetArrOfPtrs(fluxes[lev-1]),
                         refRatio(lev-1),geom[lev-1]);
      average_down_faces(GetArrOfConstPtrs(faces[lev]),
                         GetArrOfPtrs(faces[lev-1]),
                         refRatio(lev-1),geom[lev-1]);
#endif
   }

   //----------------------------------------------------------------
   // Fluxes divergence to get the velocity advection term
   for (int lev = 0; lev <= finest_level; ++lev) {

      int nGrow_divu = 4;  // Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         Real time  = getTime(lev,AmrOldTime);
         fillpatch_divu(lev,time,divu,nGrow_divu);
      }

      bool fluxes_are_area_weighted = false;
#ifdef AMREX_USE_EB
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
      //----------------------------------------------------------------
      // Use a temporary MF to hold divergence before redistribution
      int nGrow_divT = 3;
      MultiFab divTmp(grids[lev],dmap[lev],AMREX_SPACEDIM,nGrow_divT,MFInfo(),EBFactory(lev));
      divTmp.setVal(0.0);
      advFluxDivergence(lev, divTmp, 0,
                        divu,
                        GetArrOfConstPtrs(fluxes[lev]), 0,
                        GetArrOfConstPtrs(faces[lev]), 0,
                        AMREX_SPACEDIM,
                        AdvTypeVel_d.dataPtr(),
                        geom[lev], -1.0,
                        fluxes_are_area_weighted);

      divTmp.FillBoundary(geom[lev].periodicity());

      redistributeAofS(lev, m_dt,
                       divTmp, 0,
                       advData->AofS[lev], VELX,
                       ldata_p->state, VELX,
                       AMREX_SPACEDIM,
                       bcRecVel_d.dataPtr(),
                       geom[lev]);
      EB_set_covered(advData->AofS[lev],0.0);
#else
      //----------------------------------------------------------------
      // Otherwise go directly into AofS
      advFluxDivergence(lev, advData->AofS[lev], VELX,
                        divu,
                        GetArrOfConstPtrs(fluxes[lev]), 0,
                        GetArrOfConstPtrs(faces[lev]), 0,
                        AMREX_SPACEDIM,
                        AdvTypeVel_d.dataPtr(),
                        geom[lev], -1.0,
                        fluxes_are_area_weighted);
#endif
   }
}

void PeleLM::updateVelocity(std::unique_ptr<AdvanceAdvData> &advData)
{
   //----------------------------------------------------------------
   // Compute t^n divTau
   Vector<MultiFab> divtau(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), Factory(lev));
   }
   int use_density = 0;
   Real CrankNicholsonFactor = 0.5;
   computeDivTau(AmrOldTime,GetVecOfPtrs(divtau),use_density,CrankNicholsonFactor);

   //----------------------------------------------------------------
   // Get velocity forcing at half time including lagged grad P term
   int nGrow_force = 1;
   Vector<MultiFab> velForces(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      velForces[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force);
   }
   int add_gradP = 1;
   getVelForces(AmrHalfTime,GetVecOfPtrs(divtau),GetVecOfPtrs(velForces),nGrow_force,add_gradP);

   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get both old and new level_data
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

      //----------------------------------------------------------------
      // Compute provisional new velocity
      // velForce holds: 1/\rho^{n+1/2} [(gravity+...)^{n+1/2} - \nabla pi^{n} + 0.5 * divTau^{n}]
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataOld_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         auto const& vel_old  = ldataOld_p->state.const_array(mfi,VELX);
         auto const& vel_aofs = advData->AofS[lev].const_array(mfi,VELX);
         auto const& force    = velForces[lev].const_array(mfi);
         auto const& vel_new  = ldataNew_p->state.array(mfi,VELX);
         Real dt_loc = m_dt;
         amrex::ParallelFor(bx, AMREX_SPACEDIM,
         [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
            vel_new(i,j,k,n) = vel_old(i,j,k,n) + dt_loc * ( vel_aofs(i,j,k,n) + force(i,j,k,n));
         });
      }
   }
}

void PeleLM::getScalarAdvForce(std::unique_ptr<AdvanceAdvData> &advData,
                               std::unique_ptr<AdvanceDiffData> &diffData)
{
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get t^{n} data pointer
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataR_p = getLevelDataReactPtr(lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(advData->Forcing[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = ldata_p->state.const_array(mfi,DENSITY);
         auto const& rhoY    = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& T       = ldata_p->state.const_array(mfi,TEMP);
         auto const& dn      = diffData->Dn[lev].const_array(mfi,0);
         auto const& ddn     = diffData->Dn[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& r       = ldataR_p->I_R.const_array(mfi);
         auto const& extRhoY = m_extSource[lev]->const_array(mfi,FIRSTSPEC);
         auto const& extRhoH = m_extSource[lev]->const_array(mfi,RHOH);
         auto const& fY      = advData->Forcing[lev].array(mfi,0);
         auto const& fT      = advData->Forcing[lev].array(mfi,NUM_SPECIES);
         amrex::ParallelFor(bx, [rho, rhoY, T, dn, ddn, r, fY, fT, extRhoY, extRhoH, dp0dt=m_dp0dt,
                                 is_closed_ch=m_closed_chamber, do_react=m_do_react]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            buildAdvectionForcing( i, j, k, rho, rhoY, T, dn, ddn, r, extRhoY, extRhoH,
                                   dp0dt, is_closed_ch, do_react, fY, fT );
         });
      }
   }

   // Fill forcing ghost cells
   if ( advData->Forcing[0].nGrow() > 0 ) {
      fillpatch_forces(m_cur_time, GetVecOfPtrs(advData->Forcing), advData->Forcing[0].nGrow());
   }
}

void PeleLM::computeScalarAdvTerms(std::unique_ptr<AdvanceAdvData> &advData)
{

   //----------------------------------------------------------------
   // Get the BCRecs and AdvectionTypes
   auto bcRecSpec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES);
   auto bcRecSpec_d = convertToDeviceVector(bcRecSpec);
   auto AdvTypeSpec = fetchAdvTypeArray(FIRSTSPEC,NUM_SPECIES);
   auto AdvTypeSpec_d = convertToDeviceVector(AdvTypeSpec);
   auto bcRecTemp = fetchBCRecArray(TEMP,1);
   auto bcRecTemp_d = convertToDeviceVector(bcRecTemp);
   auto AdvTypeTemp = fetchAdvTypeArray(TEMP,1);
   auto AdvTypeTemp_d = convertToDeviceVector(AdvTypeTemp);
   auto bcRecRhoH = fetchBCRecArray(RHOH,1);
   auto bcRecRhoH_d = convertToDeviceVector(bcRecRhoH);
   auto AdvTypeRhoH = fetchAdvTypeArray(RHOH,1);
   auto AdvTypeRhoH_d = convertToDeviceVector(AdvTypeRhoH);

   //----------------------------------------------------------------
   // Create temporary containers
   Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         fluxes[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                  dmap[lev], NUM_SPECIES+1, 0, MFInfo(), Factory(lev));         //Species + RhoH
      }
   }

   //----------------------------------------------------------------
   // Loop over levels and get the fluxes
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get level data ptr Old
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);

      // Define edge state: Density + Species + RhoH + Temp
      int nGrow = 0;
      Array<MultiFab,AMREX_SPACEDIM> edgeState;
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         edgeState[idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                               dmap[lev], NUM_SPECIES+3, nGrow, MFInfo(), Factory(lev));
      }

#ifdef PELE_USE_EFIELD
      //----------------------------------------------------------------
      // Assemble drift and mac velocities
      ionDriftAddUmac(lev,advData);
#endif

      //----------------------------------------------------------------
      // Get divU
      int nGrow_divu = 1;  // Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         // TODO: check the number of ghost cells
         MultiFab::Copy(divu,advData->mac_divu[lev],0,0,1,nGrow_divu);
      }

      //----------------------------------------------------------------
#ifdef AMREX_USE_EB
      // Get EBFact & areafrac
      const auto& ebfact = EBFactory(lev);
      Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
      areafrac  = ebfact.getAreaFrac();
#endif

      // Get the species edge state and advection term
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& fx = fluxes[lev][0].array(mfi,0);,
                      auto const& fy = fluxes[lev][1].array(mfi,0);,
                      auto const& fz = fluxes[lev][2].array(mfi,0);)
         AMREX_D_TERM(auto const& edgex = edgeState[0].array(mfi,1);,
                      auto const& edgey = edgeState[1].array(mfi,1);,
                      auto const& edgez = edgeState[2].array(mfi,1);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& rhoY_arr  = ldata_p->state.const_array(mfi,FIRSTSPEC);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,0);

#ifdef PELE_USE_EFIELD
         // Uncharged species all at once
         bool is_velocity = false;
         bool fluxes_are_area_weighted = false;
         bool knownEdgeState = false;
         HydroUtils::ComputeFluxesOnBoxFromState(bx, NUM_SPECIES-NUM_IONS, mfi,
                                                 rhoY_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 AMREX_D_DECL(edgex,edgey,edgez), knownEdgeState,
                                                 AMREX_D_DECL(umac, vmac, wmac),
                                                 divu_arr, force_arr,
                                                 geom[lev], m_dt,
                                                 bcRecSpec, bcRecSpec_d.dataPtr(), AdvTypeSpec_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                 ebfact,
#endif
                                                 m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                 is_velocity, fluxes_are_area_weighted,
                                                 m_advection_type);

         // Ions one by one
         for ( int n = 0; n < NUM_IONS; n++) {
            auto bcRecIons = fetchBCRecArray(FIRSTSPEC+NUM_SPECIES-NUM_IONS+n,1);
            auto bcRecIons_d = convertToDeviceVector(bcRecIons);
            auto AdvTypeIons = fetchAdvTypeArray(FIRSTSPEC+NUM_SPECIES-NUM_IONS+n,1);
            auto AdvTypeIons_d = convertToDeviceVector(AdvTypeIons);
            AMREX_D_TERM(auto const& udrift = advData->uDrift[lev][0].const_array(mfi,n);,
                         auto const& vdrift = advData->uDrift[lev][1].const_array(mfi,n);,
                         auto const& wdrift = advData->uDrift[lev][2].const_array(mfi,n);)
            AMREX_D_TERM(auto const& fx_ions = fluxes[lev][0].array(mfi,NUM_SPECIES-NUM_IONS+n);,
                         auto const& fy_ions = fluxes[lev][1].array(mfi,NUM_SPECIES-NUM_IONS+n);,
                         auto const& fz_ions = fluxes[lev][2].array(mfi,NUM_SPECIES-NUM_IONS+n);)
            AMREX_D_TERM(auto const& edgex_ions = edgeState[0].array(mfi,1+NUM_SPECIES-NUM_IONS+n);,
                         auto const& edgey_ions = edgeState[1].array(mfi,1+NUM_SPECIES-NUM_IONS+n);,
                         auto const& edgez_ions = edgeState[2].array(mfi,1+NUM_SPECIES-NUM_IONS+n);)
            auto const& rhoYions_arr  = ldata_p->state.const_array(mfi,FIRSTSPEC+NUM_SPECIES-NUM_IONS+n);
            auto const& forceions_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES-NUM_IONS+n);
            HydroUtils::ComputeFluxesOnBoxFromState(bx, 1, mfi,
                                                    rhoYions_arr,
                                                    AMREX_D_DECL(fx_ions,fy_ions,fz_ions),
                                                    AMREX_D_DECL(edgex_ions,edgey_ions,edgez_ions), knownEdgeState,
                                                    AMREX_D_DECL(udrift, vdrift, wdrift),
                                                    divu_arr, forceions_arr,
                                                    geom[lev], m_dt,
                                                    bcRecIons, bcRecIons_d.dataPtr(), AdvTypeIons_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                    ebfact,
#endif
                                                    m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                    is_velocity, fluxes_are_area_weighted,
                                                    m_advection_type);
         }
#else
         bool is_velocity = false;
         bool fluxes_are_area_weighted = false;
         bool knownEdgeState = false;
         HydroUtils::ComputeFluxesOnBoxFromState(bx, NUM_SPECIES, mfi,
                                                 rhoY_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 AMREX_D_DECL(edgex,edgey,edgez), knownEdgeState,
                                                 AMREX_D_DECL(umac, vmac, wmac),
                                                 divu_arr, force_arr,
                                                 geom[lev], m_dt,
                                                 bcRecSpec, bcRecSpec_d.dataPtr(), AdvTypeSpec_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                 ebfact,
#endif
                                                 m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                 is_velocity, fluxes_are_area_weighted,
                                                 m_advection_type);
#endif
      }

      // Get edge density by summing over the species
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
         auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
#endif

         // Edge states
         for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
             const Box& ebx = amrex::surroundingNodes(bx,idim);
             auto const& rho_ed   = edgeState[idim].array(mfi,0);
             auto const& rhoY_ed  = edgeState[idim].array(mfi,1);
#ifdef AMREX_USE_EB
             if (flagfab.getType(ebx) == FabType::covered) {             // Covered boxes
                 amrex::ParallelFor(ebx, [rho_ed]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                     rho_ed(i,j,k) = 0.0;
                 });
             } else if (flagfab.getType(ebx) != FabType::regular ) {     // EB containing boxes
                 const auto& afrac = areafrac[idim]->array(mfi);
                 amrex::ParallelFor(ebx, [rho_ed, rhoY_ed, afrac]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                     rho_ed(i,j,k) = 0.0;
                     if (afrac(i,j,k) > 0.0) {                         // Uncovered faces
                         for (int n = 0; n < NUM_SPECIES; n++) {
                            rho_ed(i,j,k) += rhoY_ed(i,j,k,n);
                         }
                     }
                 });
             } else                                                     // Regular boxes
#endif
             amrex::ParallelFor(ebx, [rho_ed, rhoY_ed]
             AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                rho_ed(i,j,k) = 0.0;
                for (int n = 0; n < NUM_SPECIES; n++) {
                   rho_ed(i,j,k) += rhoY_ed(i,j,k,n);
                }
             });
         }
      }

      // Get the edge temperature
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& fx = fluxes[lev][0].array(mfi,NUM_SPECIES);,      // Put temp fluxes in place of rhoH
                      auto const& fy = fluxes[lev][1].array(mfi,NUM_SPECIES);,      // will be overwritten later
                      auto const& fz = fluxes[lev][2].array(mfi,NUM_SPECIES);)
         AMREX_D_TERM(auto const& edgex = edgeState[0].array(mfi,NUM_SPECIES+2);,
                      auto const& edgey = edgeState[1].array(mfi,NUM_SPECIES+2);,
                      auto const& edgez = edgeState[2].array(mfi,NUM_SPECIES+2);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& temp_arr  = ldata_p->state.const_array(mfi,TEMP);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES);
         bool is_velocity = false;
         bool fluxes_are_area_weighted = false;
         bool knownEdgeState = false;
         HydroUtils::ComputeFluxesOnBoxFromState(bx, 1, mfi,
                                                 temp_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 AMREX_D_DECL(edgex,edgey,edgez), knownEdgeState,
                                                 AMREX_D_DECL(umac, vmac, wmac),
                                                 divu_arr, force_arr,
                                                 geom[lev], m_dt,
                                                 bcRecTemp, bcRecTemp_d.dataPtr(), AdvTypeTemp_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                 ebfact,
#endif
                                                 m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                 is_velocity, fluxes_are_area_weighted,
                                                 m_advection_type);
      }


      // Get the edge RhoH states
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();

#ifdef AMREX_USE_EB
         auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
#endif

         for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
             const Box& ebx = amrex::surroundingNodes(bx,idim);
             auto const& rho     = edgeState[idim].const_array(mfi,0);
             auto const& rhoY    = edgeState[idim].const_array(mfi,1);
             auto const& T       = edgeState[idim].const_array(mfi,NUM_SPECIES+2);
             auto const& rhoHm   = edgeState[idim].array(mfi,NUM_SPECIES+1);
#ifdef AMREX_USE_EB
             if (flagfab.getType(ebx) == FabType::covered) {             // Covered boxes
                 amrex::ParallelFor(ebx, [rhoHm]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                     rhoHm(i,j,k) = 0.0;
                 });
             } else if (flagfab.getType(ebx) != FabType::regular ) {     // EB containing boxes
                 const auto& afrac = areafrac[idim]->array(mfi);
                 amrex::ParallelFor(ebx, [rho, rhoY, T, rhoHm, afrac]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                     if (afrac(i,j,k) <= 0.0) {                         // Covered faces
                         rhoHm(i,j,k) = 0.0;
                     } else {
                         getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
                     }
                 });
             } else                                                     // Regular boxes
#endif
             {
                 amrex::ParallelFor(ebx, [rho, rhoY, T, rhoHm]
                 AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                 {
                     getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
                 });
             }
         }
      }

      // Finally get the RhoH advection term
      // Pass the Temp forces again here, but they aren't used.
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& fx = fluxes[lev][0].array(mfi,NUM_SPECIES);,
                      auto const& fy = fluxes[lev][1].array(mfi,NUM_SPECIES);,
                      auto const& fz = fluxes[lev][2].array(mfi,NUM_SPECIES);)
         AMREX_D_TERM(auto const& edgex = edgeState[0].array(mfi,NUM_SPECIES+1);,
                      auto const& edgey = edgeState[1].array(mfi,NUM_SPECIES+1);,
                      auto const& edgez = edgeState[2].array(mfi,NUM_SPECIES+1);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& rhoh_arr  = ldata_p->state.const_array(mfi,RHOH);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES);
         bool is_velocity = false;
         bool fluxes_are_area_weighted = false;
         bool knownEdgeState = true;
         HydroUtils::ComputeFluxesOnBoxFromState(bx, 1, mfi,
                                                 rhoh_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 AMREX_D_DECL(edgex,edgey,edgez), knownEdgeState,
                                                 AMREX_D_DECL(umac, vmac, wmac),
                                                 divu_arr, force_arr,
                                                 geom[lev], m_dt,
                                                 bcRecRhoH, bcRecRhoH_d.dataPtr(), AdvTypeRhoH_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                 ebfact,
#endif
                                                 m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                 is_velocity, fluxes_are_area_weighted,
                                                 m_advection_type);
      }
#ifdef AMREX_USE_EB
      EB_set_covered_faces(GetArrOfPtrs(fluxes[lev]),0.);
#endif
   }

   //----------------------------------------------------------------
   // Average down fluxes to ensure C/F consistency
   for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
      EB_average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                            GetArrOfPtrs(fluxes[lev-1]),
                            refRatio(lev-1),geom[lev-1]);
#else
      average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                         GetArrOfPtrs(fluxes[lev-1]),
                         refRatio(lev-1),geom[lev-1]);
#endif
   }

   //----------------------------------------------------------------
   // If balances are required, compute face domain integrals
   // using level 0 since we've averaged down the fluxes already
   if (m_sdcIter == m_nSDCmax) {
      if (m_do_massBalance) addMassFluxes(GetArrOfConstPtrs(fluxes[0]),geom[0]);
      if (m_do_energyBalance) addRhoHFluxes(GetArrOfConstPtrs(fluxes[0]),geom[0]);
      if (m_do_speciesBalance) addRhoYFluxes(GetArrOfConstPtrs(fluxes[0]),geom[0]);
   }
   // Compute face domain integral for U at every SDC iteration
   addUmacFluxes(advData, geom[0]);

#ifdef PELE_USE_EFIELD
   if (m_do_extraEFdiags) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         for (int n = 0; n < NUM_IONS; ++n) {
            int spec_idx = NUM_SPECIES - NUM_IONS + n;
            Array<std::unique_ptr<MultiFab>,AMREX_SPACEDIM> ionFlux;
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++ ) {
               ionFlux[idim].reset(new MultiFab(fluxes[lev][idim],amrex::make_alias,spec_idx,1));
            }
            average_face_to_cellcenter(*m_ionsFluxes[lev],n*AMREX_SPACEDIM,
                                       GetArrOfConstPtrs(ionFlux));
         }
      }
   }
#endif

   //----------------------------------------------------------------
   // Fluxes divergence to get the scalars advection term
   auto AdvTypeAll = fetchAdvTypeArray(FIRSTSPEC,NUM_SPECIES+1); // Species+RhoH
   auto AdvTypeAll_d = convertToDeviceVector(AdvTypeAll);
   for (int lev = 0; lev <= finest_level; ++lev) {

      int nGrow_divu = 1;  // TODO EB Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         Real time  = getTime(lev,AmrOldTime);
         fillpatch_divu(lev,time,divu,nGrow_divu);
      }

      bool fluxes_are_area_weighted = false;
#ifdef AMREX_USE_EB
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
      //----------------------------------------------------------------
      // Use a temporary MF to hold divergence before redistribution
      int nGrow_divTmp= 3;
      MultiFab divTmp(grids[lev],dmap[lev],NUM_SPECIES+1,nGrow_divTmp,MFInfo(),EBFactory(lev));
      divTmp.setVal(0.0);
      advFluxDivergence(lev, divTmp, 0,
                        divu,
                        GetArrOfConstPtrs(fluxes[lev]), 0,
                        GetArrOfConstPtrs(fluxes[lev]), 0, // This will not be used since none of rhoY/rhoH in convective
                        NUM_SPECIES+1,
                        AdvTypeAll_d.dataPtr(),
                        geom[lev], -1.0,
                        fluxes_are_area_weighted);

      divTmp.FillBoundary(geom[lev].periodicity());

      // Need separate redistribution of rhoYs/rhoH
      redistributeAofS(lev, m_dt,
                       divTmp, 0,
                       advData->AofS[lev], FIRSTSPEC,
                       ldata_p->state, FIRSTSPEC,
                       NUM_SPECIES,
                       bcRecSpec_d.dataPtr(),
                       geom[lev]);

      redistributeAofS(lev, m_dt,
                       divTmp, NUM_SPECIES,
                       advData->AofS[lev], RHOH,
                       ldata_p->state, RHOH,
                       1,
                       bcRecRhoH_d.dataPtr(),
                       geom[lev]);
      EB_set_covered(advData->AofS[lev],0.0);
#else
      //----------------------------------------------------------------
      // Otherwise go directly into AofS
      advFluxDivergence(lev, advData->AofS[lev], FIRSTSPEC,
                        divu,
                        GetArrOfConstPtrs(fluxes[lev]), 0,
                        GetArrOfConstPtrs(fluxes[lev]), 0, // This will not be used since none of rhoY/rhoH in convective
                        NUM_SPECIES+1,
                        AdvTypeAll_d.dataPtr(),
                        geom[lev], -1.0,
                        fluxes_are_area_weighted);
#endif
   }

   //----------------------------------------------------------------
   // Sum over the species AofS to get the density advection term
   for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(advData->AofS[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& bx = mfi.tilebox();
         auto const& aofrho = advData->AofS[lev].array(mfi,DENSITY);
         auto const& aofrhoY = advData->AofS[lev].const_array(mfi,FIRSTSPEC);
         amrex::ParallelFor(bx, [aofrho, aofrhoY]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            aofrho(i,j,k) = 0.0;
            for (int n = 0; n < NUM_SPECIES; n++) {
               aofrho(i,j,k) += aofrhoY(i,j,k,n);
            }
         });
      }
   }
}

void PeleLM::updateDensity(std::unique_ptr<AdvanceAdvData> &advData)
{
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get level data ptr
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNew_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& bx = mfi.tilebox();
         auto const& rhoOld_arr  = ldataOld_p->state.const_array(mfi,DENSITY);
         auto const& rhoNew_arr  = ldataNew_p->state.array(mfi,DENSITY);
         auto const& a_of_rho    = advData->AofS[lev].const_array(mfi,DENSITY);
         auto const& ext_rho     = m_extSource[lev]->const_array(mfi,DENSITY);
         amrex::ParallelFor(bx, [rhoOld_arr, rhoNew_arr, a_of_rho, ext_rho, dt=m_dt]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            rhoNew_arr(i,j,k) = rhoOld_arr(i,j,k) + dt * (a_of_rho(i,j,k) + ext_rho(i,j,k));
         });
      }
   }
   averageDownDensity(AmrNewTime);
}

void PeleLM::computePassiveAdvTerms(std::unique_ptr<AdvanceAdvData> &advData,
                                    int state_comp,
                                    int ncomp)
{
   //----------------------------------------------------------------
   // Get the BCRecs and AdvectionTypes
   auto bcRecPass = fetchBCRecArray(state_comp,ncomp);
   auto bcRecPass_d = convertToDeviceVector(bcRecPass);
   auto AdvTypePass = fetchAdvTypeArray(state_comp,ncomp);
   auto AdvTypePass_d = convertToDeviceVector(AdvTypePass);

   //----------------------------------------------------------------
   // Create temporary containers
   Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
   Vector<Array<MultiFab,AMREX_SPACEDIM> > edgeState(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         fluxes[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                  dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
         edgeState[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
      }
   }

   //----------------------------------------------------------------
   // Loop over levels and get the fluxes
   for (int lev = 0; lev <= finest_level; ++lev) {
      // Get level data ptr Old
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);

      //----------------------------------------------------------------
      // Get divU
      int nGrow_divu = 1;  // Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         // TODO: check the number of ghost cells
         MultiFab::Copy(divu,advData->mac_divu[lev],0,0,1,nGrow_divu);
      }

      //----------------------------------------------------------------
#ifdef AMREX_USE_EB
      // Get EBFact & areafrac
      const auto& ebfact = EBFactory(lev);
      Array< const MultiCutFab*,AMREX_SPACEDIM> areafrac;
      areafrac  = ebfact.getAreaFrac();
#endif

      // Get the passive variables edge state and advection term
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& fx = fluxes[lev][0].array(mfi,0);,
                      auto const& fy = fluxes[lev][1].array(mfi,0);,
                      auto const& fz = fluxes[lev][2].array(mfi,0);)
         AMREX_D_TERM(auto const& edgex = edgeState[lev][0].array(mfi,0);,
                      auto const& edgey = edgeState[lev][1].array(mfi,0);,
                      auto const& edgez = edgeState[lev][2].array(mfi,0);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& pass_arr  = ldata_p->state.const_array(mfi,state_comp);
         // TODO: Find way to include diffusive forces for passive scalars that diffuse
         auto const& force_arr = m_extSource[lev]->const_array(mfi,state_comp);
         bool is_velocity = false;
         bool fluxes_are_area_weighted = false;
         bool knownEdgeState = false;
         HydroUtils::ComputeFluxesOnBoxFromState(bx, ncomp, mfi,
                                                 pass_arr,
                                                 AMREX_D_DECL(fx,fy,fz),
                                                 AMREX_D_DECL(edgex,edgey,edgez), knownEdgeState,
                                                 AMREX_D_DECL(umac, vmac, wmac),
                                                 divu_arr, force_arr,
                                                 geom[lev], m_dt,
                                                 bcRecPass, bcRecPass_d.dataPtr(), AdvTypePass_d.dataPtr(),
#ifdef AMREX_USE_EB
                                                 ebfact,
#endif
                                                 m_Godunov_ppm, m_Godunov_ForceInTrans,
                                                 is_velocity, fluxes_are_area_weighted,
                                                 m_advection_type);
      }
#ifdef AMREX_USE_EB
      EB_set_covered_faces(GetArrOfPtrs(fluxes[lev]),0.);
#endif
   }
   //----------------------------------------------------------------
   // Average down fluxes and edge state to ensure C/F consistency
   for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
      EB_average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                            GetArrOfPtrs(fluxes[lev-1]),
                            refRatio(lev-1),geom[lev-1]);
      EB_average_down_faces(GetArrOfConstPtrs(edgeState[lev]),
                            GetArrOfPtrs(edgeState[lev-1]),
                            refRatio(lev-1),geom[lev-1]);
#else
      average_down_faces(GetArrOfConstPtrs(fluxes[lev]),
                         GetArrOfPtrs(fluxes[lev-1]),
                         refRatio(lev-1),geom[lev-1]);
      average_down_faces(GetArrOfConstPtrs(edgeState[lev]),
                         GetArrOfPtrs(edgeState[lev-1]),
                         refRatio(lev-1),geom[lev-1]);
#endif
   }
   //----------------------------------------------------------------
   // Fluxes divergence to get the scalars advection term
   auto AdvTypeAll = fetchAdvTypeArray(state_comp,ncomp);
   auto AdvTypeAll_d = convertToDeviceVector(AdvTypeAll);
   for (int lev = 0; lev <= finest_level; ++lev) {

      int nGrow_divu = 1;  // TODO EB Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         Real time  = getTime(lev,AmrOldTime);
         fillpatch_divu(lev,time,divu,nGrow_divu);
      }

      bool fluxes_are_area_weighted = false;
#ifdef AMREX_USE_EB
      auto ldata_p = getLevelDataPtr(lev,AmrOldTime);
      //----------------------------------------------------------------
      // Use a temporary MF to hold divergence before redistribution
      int nGrow_divTmp= 3;
      MultiFab divTmp(grids[lev],dmap[lev],ncomp,nGrow_divTmp,MFInfo(),EBFactory(lev));
      divTmp.setVal(0.0);
      advFluxDivergence(lev, divTmp, 0,
                        divu,
                        GetArrOfConstPtrs(fluxes[lev]), 0,
                        GetArrOfConstPtrs(edgeState[lev]), 0,
                        ncomp,
                        AdvTypeAll_d.dataPtr(),
                        geom[lev], -1.0,
                        fluxes_are_area_weighted);

      divTmp.FillBoundary(geom[lev].periodicity());

      redistributeAofS(lev, m_dt,
                       divTmp, 0,
                       advData->AofS[lev], state_comp,
                       ldata_p->state, state_comp,
                       ncomp,
                       bcRecPass_d.dataPtr(),
                       geom[lev]);
#else
      //----------------------------------------------------------------
      // Otherwise go directly into AofS
      advFluxDivergence(lev, advData->AofS[lev], state_comp,
                        divu,
                        GetArrOfConstPtrs(fluxes[lev]), 0,
                        GetArrOfConstPtrs(edgeState[lev]), 0,
                        ncomp,
                        AdvTypeAll_d.dataPtr(),
                        geom[lev], -1.0,
                        fluxes_are_area_weighted);
#endif
   }
   // TODO: This assumes passive variables have no diffusive fluxes
   updateScalarComp(advData, state_comp, ncomp);
}

void PeleLM::updateScalarComp(std::unique_ptr<AdvanceAdvData> &advData,
                              int state_comp,
                              int ncomp)
{
   for (int lev = 0; lev <= finest_level; ++lev) {

      // Get level data ptr
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNew_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& bx = mfi.tilebox();
         auto const& old_arr  = ldataOld_p->state.const_array(mfi,state_comp);
         auto const& new_arr  = ldataNew_p->state.array(mfi,state_comp);
         auto const& a_of_s    = advData->AofS[lev].const_array(mfi,state_comp);
         auto const& ext     = m_extSource[lev]->const_array(mfi,state_comp);
         amrex::ParallelFor(bx, ncomp, [old_arr, new_arr, a_of_s, ext, dt=m_dt]
         AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
         {
           new_arr(i,j,k,n) = old_arr(i,j,k,n) + dt * (a_of_s(i,j,k,n) + ext(i,j,k,n));
         });
      }
   }
   averageDown(AmrNewTime, state_comp, ncomp);
}
