#include <PeleLM.H>
#include <PeleLMUtils.H>
#include <Godunov.H>

using namespace amrex;

void PeleLM::Advance(int is_initIter) {
   BL_PROFILE_VAR("PeleLM::Advance()", Advance);

   // Start timing current time step
   Real strt_time = ParallelDescriptor::second();

   //----------------------------------------------------------------
   // Copy old <- new state
   copyStateNewToOld();
   copyPressNewToOld();
   //----------------------------------------------------------------

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
   advData.reset(new AdvanceAdvData(finest_level, grids, dmap, m_factory, m_incompressible, 
                                    m_nGrowAdv, m_nGrowMAC));
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
      fillPatchState(AmrNewTime);
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
   //----------------------------------------------------------------
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
   //----------------------------------------------------------------
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
   // Get advection SDC forcing
   getScalarAdvForce(advData);

   // Get AofS (\nabla \cdot (\rho Y Umac))^{n+1/2}
   // and density = \sum_k AofS_k
   computeScalarAdvTerms(advData);
   
   // Compute \rho^{np1,k} and fillpatch new density
   updateDensity(advData);
   fillPatchDensity(AmrNewTime);
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar diffusion
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Reaction
   //----------------------------------------------------------------
   // TODO: right now just add Adv and Diff
   for (int lev = 0; lev <= finest_level; ++lev) { 

      // Get level data ptr
      auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
      auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldataNew_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& bx = mfi.tilebox();
         auto const& rhoYo_arr  = ldataOld_p->species.const_array(mfi);
         auto const& rhoYn_arr  = ldataNew_p->species.array(mfi);
         auto const& rhoHo_arr  = ldataOld_p->rhoh.const_array(mfi);
         auto const& rhoHn_arr  = ldataNew_p->rhoh.array(mfi);
         auto const& a_of_rhoY   = advData->AofS[lev].const_array(mfi,FIRSTSPEC);
         auto const& a_of_rhoH   = advData->AofS[lev].const_array(mfi,RHOH);
         amrex::ParallelFor(bx, [=,dt = m_dt]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            for (int n = 0; n < NUM_SPECIES; n++ ) {
               rhoYn_arr(i,j,k,n) = rhoYo_arr(i,j,k,n) + dt * a_of_rhoY(i,j,k,n);
            }
            rhoHn_arr(i,j,k) = rhoHo_arr(i,j,k) + dt * a_of_rhoH(i,j,k);
         });
      }
   }

   //----------------------------------------------------------------
   // Wrap it up
   //----------------------------------------------------------------
   setTemperature(AmrNewTime);
}

void PeleLM::updateVelocityDiffusion(std::unique_ptr<AdvanceDiffData> &diffData)
{
}
