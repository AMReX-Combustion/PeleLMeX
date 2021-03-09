#include <PeleLM.H>
#include <PeleLM_K.H>
#include <PeleLMUtils.H>
#include <Godunov.H>

using namespace amrex;

void PeleLM::computeVelocityAdvTerm(std::unique_ptr<AdvanceAdvData> &advData)
{

   //----------------------------------------------------------------
   // Get viscous forces
   int nGrow_force = 1;
   Vector<MultiFab> divtau(finest_level+1);
   Vector<MultiFab> velForces(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      divtau[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, 0);
      velForces[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force);
   }
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
      // Compute the velocity advective term

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& vel_arr   = ldata_p->velocity.const_array(mfi);
         auto const& force_arr = velForces[lev].const_array(mfi);
         auto const& aofs_vel  = advData->AofS[lev].array(mfi,VELX);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;   
         FArrayBox tmpfab(amrex::grow(bx,1), AMREX_SPACEDIM*n_tmp_fac+1);

#ifdef AMREX_USE_EB
#else
         godunov::compute_godunov_advection(bx, AMREX_SPACEDIM,
                                            aofs_vel, vel_arr,
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

void PeleLM::getScalarAdvForce(std::unique_ptr<AdvanceAdvData> &advData)
{
   // TODO: forcing for advection
   for (int lev = 0; lev <= finest_level; ++lev) {   
      advData->Forcing[lev].setVal(0.0);
   }
}

void PeleLM::computeScalarAdvTerms(std::unique_ptr<AdvanceAdvData> &advData)
{

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
      int nGrow_divu = 4;  // Why incflo use 4 ?
      MultiFab divu(grids[lev],dmap[lev],1,nGrow_divu,MFInfo(),Factory(lev));
      if (m_incompressible) {
         divu.setVal(0.0);
      } else {
         // TODO: this should be mac_divu
         Real time  = getTime(lev,AmrOldTime);
         fillpatch_divu(lev,time,divu,nGrow_divu);
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
         AMREX_D_TERM(auto const& edgex = edgeState[0].array(mfi,1);,
                      auto const& edgey = edgeState[1].array(mfi,1);,
                      auto const& edgez = edgeState[2].array(mfi,1);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& rhoY_arr  = ldata_p->species.const_array(mfi);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,0);
         auto const& aofs      = advData->AofS[lev].array(mfi,FIRSTSPEC);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;   
         FArrayBox tmpfab(amrex::grow(bx,1), NUM_SPECIES*n_tmp_fac+1);

#ifdef AMREX_USE_EB
#else
         godunov::compute_godunov_advection(bx, NUM_SPECIES,
                                            aofs, rhoY_arr,
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

      // Get the density AofS and edge density by summing over the species
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();

         // Advection term
         auto const& adv_rho   = advData->AofS[lev].array(mfi,DENSITY);
         auto const& adv_rhoY  = advData->AofS[lev].array(mfi,FIRSTSPEC);
         amrex::ParallelFor(bx, [ adv_rho, adv_rhoY ]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            adv_rho(i,j,k) = 0.0;
            for (int n = 0; n < NUM_SPECIES; n++) {
               adv_rho(i,j,k) += adv_rhoY(i,j,k,n); 
            }
         });

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
      // Discard the AofS, not used
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& edgex = edgeState[0].array(mfi,NUM_SPECIES+2);,
                      auto const& edgey = edgeState[1].array(mfi,NUM_SPECIES+2);,
                      auto const& edgez = edgeState[2].array(mfi,NUM_SPECIES+2);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& temp_arr  = ldata_p->temp.const_array(mfi);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES-1);
         auto const& aofs      = advData->AofS[lev].array(mfi,TEMP);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;   
         FArrayBox tmpfab(amrex::grow(bx,1), 1*n_tmp_fac+1);

#ifdef AMREX_USE_EB
#else
         godunov::compute_godunov_advection(bx, 1,
                                            aofs, temp_arr,
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
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

         Box const& bx = mfi.tilebox();
         AMREX_D_TERM(auto const& umac = advData->umac[lev][0].const_array(mfi);,
                      auto const& vmac = advData->umac[lev][1].const_array(mfi);,
                      auto const& wmac = advData->umac[lev][2].const_array(mfi);)
         AMREX_D_TERM(auto const& edgex = edgeState[0].array(mfi,NUM_SPECIES+1);,
                      auto const& edgey = edgeState[1].array(mfi,NUM_SPECIES+1);,
                      auto const& edgez = edgeState[2].array(mfi,NUM_SPECIES+1);)
         auto const& divu_arr  = divu.const_array(mfi);
         auto const& rhoh_arr  = ldata_p->rhoh.const_array(mfi);
         auto const& force_arr = advData->Forcing[lev].const_array(mfi,NUM_SPECIES);
         auto const& aofs      = advData->AofS[lev].array(mfi,RHOH);
         // Temporary fab used in Godunov. TODO ngrow should be > 1 for EB.
         int n_tmp_fac = (AMREX_SPACEDIM == 2) ? 10 : 14;   
         FArrayBox tmpfab(amrex::grow(bx,1), 1*n_tmp_fac+1);

#ifdef AMREX_USE_EB
#else
         godunov::compute_godunov_advection(bx, 1,
                                            aofs, rhoh_arr,
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
