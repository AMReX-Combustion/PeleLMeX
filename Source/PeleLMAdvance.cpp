#include <PeleLM.H>
#include <PeleLMUtils.H>
#include <Godunov.H>

using namespace amrex;

void PeleLM::Advance(int is_initIter) {
   BL_PROFILE_VAR("PeleLM::Advance()", Advance);

   // Start timing current time step
   Real strt_time = ParallelDescriptor::second();

   //----------------------------------------------------------------
   // TIME
   // Compute time-step size
   m_dt = computeDt(is_initIter);

   // Update time vectors
   for(int lev = 0; lev <= finest_level; lev++)
   {   
       m_t_old[lev] = m_cur_time;
       m_t_new[lev] = m_cur_time + m_dt;
   } 
   //----------------------------------------------------------------

   if ( m_verbose ) {
      amrex::Print() << " STEP [" << m_nstep << "] - Time: " << m_cur_time << ", dt " << m_dt << "\n";
   }

   //----------------------------------------------------------------
   // Data for the advance
   std::unique_ptr<AdvanceDiffData> diffData;
   diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory, m_nGrowAdv, m_use_wbar));
   std::unique_ptr<AdvanceAdvData> advData;
   advData.reset(new AdvanceAdvData(finest_level, grids, dmap, m_factory, m_incompressible, m_nGrowMAC));
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Advance setup
   // fillpatch the t^{n} data
   fillPatchState(AmrOldTime);
   
   // compute t^n data 
   calcViscosity(AmrOldTime);
   if (! m_incompressible ) {
      calcDiffusivity(AmrOldTime);
   }
   // TODO : typical values
   // TODO : check dt
   // TODO : floor_species

   if (! m_incompressible ) {
      setThermoPress(AmrOldTime);
      computeDifferentialDiffusionTerms(AmrOldTime,diffData->Dn);
   }

   copyTransportOldToNew();
   if (! m_incompressible ) {
      copyDiffusionOldToNew(diffData);
   }

   // TODO : handle reaction ghost cells
   //----------------------------------------------------------------
   
   //----------------------------------------------------------------
   // Scalar advance
   if ( m_incompressible ) {

      // Still need to get face velocities
      predictVelocity(advData,diffData);
      // MAC-project face velocities
      macProject(AmrOldTime,advData,{});
   } else {

      // SDC iterations
      for (int sdc_iter = 1; sdc_iter <= m_nSDCmax; ++sdc_iter ) {
         oneSDC(sdc_iter,advData,diffData);
      }
      // Post SDC
      calcDiffusivity(AmrNewTime);
      setThermoPress(AmrNewTime);
      if (m_has_divu) {
         int is_initialization = 0;
         int do_avgDown = 1;
         calcDivU(is_initialization,do_avgDown,AmrNewTime);
      }
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Velocity advance
   // Re-evaluate viscosity only if scalar updated
   if (!m_incompressible) calcViscosity(AmrNewTime);

   // Compute t^{n+1/2} velocity advection term
   computeVelocityAdvTerm(advData);
   
   // Compute provisional new velocity for diffusion solve RHS
   // U^{np1**} = U^{n} - dt*AofS^{n+1/2} - dt/rho^{n+1/2} \nabla \pi^{n-1/2} + dt/rho^{n+1/2} * F^{n+1/2}
   updateVelocity(is_initIter,advData);

   // Semi-implicit CN diffusion solve to get U^{np1*}
   // updateVelocityDiffusion(diffData);

   // Nodal projection to get constrained U^{np1} and new pressure \pi^{n+1/2}
   const TimeStamp rhoTime = AmrHalfTime;
   velocityProjection(is_initIter,rhoTime,m_dt);
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Wrapup advance
   // Timing current time step
   if (m_verbose > 0)
   {
      Real run_time = ParallelDescriptor::second() - strt_time;
      ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << " >> PeleLM::Advance() " << run_time << "\n";
   }
}

void PeleLM::oneSDC(int sdcIter,
                    std::unique_ptr<AdvanceAdvData> &advData,
                    std::unique_ptr<AdvanceDiffData> &diffData)
{

   //----------------------------------------------------------------
   // Update t^{n+1,k} transport/Dnp1/divU

   // For the first SDC, we already copied old -> new
   if (sdcIter > 1) {
      calcDiffusivity(AmrNewTime);
      computeDifferentialDiffusionTerms(AmrNewTime,diffData->Dnp1);
      int is_initialization = 0;
      int do_avgDown = 1;
      calcDivU(is_initialization,do_avgDown,AmrNewTime);
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Get u MAC
   
   // Predict face velocity with Godunov
   predictVelocity(advData,diffData);

   // Create S^{n+1/2} by averaging old and new

   // Add chi_increment

   // MAC projection
   macProject(AmrOldTime,advData,{});
   //----------------------------------------------------------------
   
   //----------------------------------------------------------------
   // Scalar advections
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar diffusion
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Reaction
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Wrap it up
   //----------------------------------------------------------------
}

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
         computeDivergence(divu,GetArrOfConstPtrs(advData->umac[lev]),geom[lev]);
         divu.FillBoundary(geom[lev].periodicity());
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

void PeleLM::updateVelocityDiffusion(std::unique_ptr<AdvanceDiffData> &diffData)
{
}
