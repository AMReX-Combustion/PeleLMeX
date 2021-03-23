#include <PeleLM.H>
#include <PeleLM_K.H>
#include <PeleLMUtils.H>
#include <Godunov.H>

using namespace amrex;

void PeleLM::computeVelocityAdvTerm(std::unique_ptr<AdvanceAdvData> &advData)
{

   //----------------------------------------------------------------
   // Create temporary containers
   int nGrow_force = 1;
   Vector<MultiFab> divtau(finest_level+1);
   Vector<MultiFab> velForces(finest_level+1);
   Vector<Array<MultiFab,AMREX_SPACEDIM> > fluxes(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0, MFInfo(), Factory(lev));
      velForces[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force, MFInfo(), Factory(lev));
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         fluxes[lev][idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
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
      // Compute the velocity fluxes
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& fx = fluxes[lev][0].array(mfi);,
                      auto const& fy = fluxes[lev][1].array(mfi);,
                      auto const& fz = fluxes[lev][2].array(mfi);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& vel_arr   = ldata_p->velocity.const_array(mfi);
         auto const& force_arr = velForces[lev].const_array(mfi);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;   
         FArrayBox tmpfab(amrex::grow(bx,1), AMREX_SPACEDIM*n_tmp_fac+1);

#ifdef AMREX_USE_EB
         // TODO
#else
         godunov::compute_godunov_fluxes(bx, 0, AMREX_SPACEDIM,
                                         AMREX_D_DECL(fx,fy,fz), vel_arr,
                                         AMREX_D_DECL(umac, vmac, wmac),
                                         force_arr, divu_arr, m_dt,
                                         bcRecVel_d.dataPtr(),
                                         AdvTypeVel_d.dataPtr(),
                                         tmpfab.dataPtr(),m_Godunov_ppm,
                                         m_Godunov_ForceInTrans,
                                         geom[lev], true);
#endif
      }
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
   // Fluxes divergence to get the velocity advection term
   advFluxDivergence(advData->AofS,VELX,GetVecOfArrOfPtrs(fluxes),0,
                     GetVecOfArrOfPtrs(advData->umac),AMREX_SPACEDIM,AdvTypeVel_d.dataPtr(),-1.0);
}

void PeleLM::updateVelocity(int is_init,
                            std::unique_ptr<AdvanceAdvData> &advData)
{

   //----------------------------------------------------------------
   // Compute t^n divTau
   Vector<MultiFab> divtau(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataOld_p->velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         auto const& vel_old  = ldataOld_p->velocity.const_array(mfi);
         auto const& vel_aofs = advData->AofS[lev].const_array(mfi,VELX);
         auto const& force    = velForces[lev].const_array(mfi);
         auto const& vel_new  = ldataNew_p->velocity.array(mfi);
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

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(advData->Forcing[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto const& rho     = ldata_p->density.const_array(mfi);
         auto const& rhoY    = ldata_p->species.const_array(mfi);
         auto const& T       = ldata_p->temp.const_array(mfi);
         auto const& dn      = diffData->Dn[lev].const_array(mfi,0);
         auto const& ddn     = diffData->Dn[lev].const_array(mfi,NUM_SPECIES+1);
         auto const& r       = ldata_p->temp.const_array(mfi); // TODO get_new_data(RhoYdot_Type).array(mfi,0);
         auto const& fY      = advData->Forcing[lev].array(mfi,0);
         auto const& fT      = advData->Forcing[lev].array(mfi,NUM_SPECIES);
         Real        dp0dt_d = 0.0; // TODO dp0dt;
         int     closed_ch_d = 0;   // TODO closed_chamber;

         amrex::ParallelFor(bx, [rho, rhoY, T, dn, ddn, r, fY, fT, dp0dt_d, closed_ch_d]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            buildAdvectionForcing( i, j, k, rho, rhoY, T, dn, ddn, r, dp0dt_d, closed_ch_d, fY, fT );
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
      int nGrow = 0; // TODO more needed for EB ?
      Array<MultiFab,AMREX_SPACEDIM> edgeState;
      for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
         edgeState[idim].define(amrex::convert(grids[lev],IntVect::TheDimensionVector(idim)),
                                               dmap[lev], NUM_SPECIES+3, nGrow, MFInfo(), Factory(lev));
      }

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

      // Get the species edge state and advection term
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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
         auto const& rhoY_arr  = ldata_p->species.const_array(mfi);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,0);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;
         FArrayBox tmpfab(amrex::grow(bx,1), NUM_SPECIES*n_tmp_fac+1);

#ifdef AMREX_USE_EB
         //TODO
#else
         godunov::compute_godunov_fluxes(bx, 0, NUM_SPECIES,
                                         AMREX_D_DECL(fx,fy,fz), rhoY_arr,
                                         AMREX_D_DECL(umac, vmac, wmac),
                                         force_arr, divu_arr,
                                         AMREX_D_DECL(edgex, edgey, edgez), false,
                                         m_dt,
                                         bcRecSpec_d.dataPtr(),
                                         AdvTypeSpec_d.dataPtr(),
                                         tmpfab.dataPtr(),m_Godunov_ppm,
                                         m_Godunov_ForceInTrans,
                                         geom[lev]);
#endif
      }

      // Get edge density by summing over the species
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();

         // Edge states
         for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
            const Box& ebx = amrex::surroundingNodes(bx,idim);
            auto const& rho_ed   = edgeState[idim].array(mfi,0);
            auto const& rhoY_ed  = edgeState[idim].array(mfi,1);
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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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
         auto const& temp_arr  = ldata_p->temp.const_array(mfi);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;
         FArrayBox tmpfab(amrex::grow(bx,1), 1*n_tmp_fac+1);

#ifdef AMREX_USE_EB
         //TODO
#else
         godunov::compute_godunov_fluxes(bx, 0, 1,
                                         AMREX_D_DECL(fx,fy,fz), temp_arr,
                                         AMREX_D_DECL(umac, vmac, wmac),
                                         force_arr, divu_arr,
                                         AMREX_D_DECL(edgex, edgey, edgez), false,
                                         m_dt,
                                         bcRecTemp_d.dataPtr(),
                                         AdvTypeTemp_d.dataPtr(),
                                         tmpfab.dataPtr(),m_Godunov_ppm,
                                         m_Godunov_ForceInTrans,
                                         geom[lev]);
#endif
      }

      // Get the edge RhoH states
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();

         for (int idim = 0; idim <AMREX_SPACEDIM; idim++) {
            const Box& ebx = amrex::surroundingNodes(bx,idim);
            auto const& rho     = edgeState[idim].const_array(mfi,0);
            auto const& rhoY    = edgeState[idim].const_array(mfi,1);
            auto const& T       = edgeState[idim].const_array(mfi,NUM_SPECIES+2);
            auto const& rhoHm   = edgeState[idim].array(mfi,NUM_SPECIES+1);
            amrex::ParallelFor(ebx, [rho, rhoY, T, rhoHm]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                getRHmixGivenTY( i, j, k, rho, rhoY, T, rhoHm );
            });
         }
      }

      // Finally get the RhoH advection term
      // Pass the Temp forces again here, but they aren't used.
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

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
         auto const& rhoh_arr  = ldata_p->rhoh.const_array(mfi);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;
         FArrayBox tmpfab(amrex::grow(bx,1), 1*n_tmp_fac+1);

#ifdef AMREX_USE_EB
         //TODO
#else
         godunov::compute_godunov_fluxes(bx, 0, 1,
                                         AMREX_D_DECL(fx,fy,fz), rhoh_arr,
                                         AMREX_D_DECL(umac, vmac, wmac),
                                         force_arr, divu_arr,
                                         AMREX_D_DECL(edgex, edgey, edgez), true,
                                         m_dt,
                                         bcRecRhoH_d.dataPtr(),
                                         AdvTypeRhoH_d.dataPtr(),
                                         tmpfab.dataPtr(),m_Godunov_ppm,
                                         m_Godunov_ForceInTrans,
                                         geom[lev]);
#endif
      }
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
   // Fluxes divergence to get the scalars advection term
   auto AdvTypeAll = fetchAdvTypeArray(FIRSTSPEC,NUM_SPECIES+1); // Species+RhoH
   auto AdvTypeAll_d = convertToDeviceVector(AdvTypeAll);
   advFluxDivergence(advData->AofS,FIRSTSPEC,GetVecOfArrOfPtrs(fluxes),0,
                     GetVecOfArrOfPtrs(advData->umac),NUM_SPECIES+1,AdvTypeAll_d.dataPtr(),-1.0);

   //----------------------------------------------------------------
   // Sum over the species AofS to get the density advection term
   for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef _OPENMP
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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNew_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& bx = mfi.tilebox();
         auto const& rhoOld_arr  = ldataOld_p->density.const_array(mfi);
         auto const& rhoNew_arr  = ldataNew_p->density.array(mfi);
         auto const& a_of_rho    = advData->AofS[lev].const_array(mfi,DENSITY);
         amrex::ParallelFor(bx, [rhoOld_arr, rhoNew_arr, a_of_rho,dt=m_dt]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            rhoNew_arr(i,j,k) = rhoOld_arr(i,j,k) + dt * a_of_rho(i,j,k);
         });
      }
   }
}
