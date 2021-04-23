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
   m_dt = computeDt(is_initIter,AmrOldTime);

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
   // Data for the advance, only live for the duration of the advance
   std::unique_ptr<AdvanceDiffData> diffData;
   diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory, m_nGrowAdv, m_use_wbar));
   std::unique_ptr<AdvanceAdvData> advData;
   advData.reset(new AdvanceAdvData(finest_level, grids, dmap, m_factory, m_incompressible,
                                    m_nGrowAdv, m_nGrowMAC));
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Advance setup
   // fillpatch the t^{n} data
   averageDownState(AmrOldTime);
   fillPatchState(AmrOldTime);

   // compute t^{n} data
   calcViscosity(AmrOldTime);
   if (! m_incompressible ) {
      calcDiffusivity(AmrOldTime);
#ifdef PLM_USE_EFIELD
      poissonSolveEF(AmrOldTime);
#endif
   }
   // TODO : typical values
   // TODO : check dt
   // TODO : floor_species

   if (! m_incompressible ) {
      floorSpecies(AmrOldTime);
      setThermoPress(AmrOldTime);
      computeDifferentialDiffusionTerms(AmrOldTime,diffData);
   }

   // Initialize terms t^{np1,k} from t^{n}
   copyTransportOldToNew();
   if (! m_incompressible ) {
      copyDiffusionOldToNew(diffData);
#ifdef PLM_USE_EFIELD
      ionDriftVelocity(advData);
#endif
   }

   // TODO : handle reaction ghost cells
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar advance
   if ( m_incompressible ) {

      // Still need to get face velocities ...
      predictVelocity(advData,diffData);

      // ... and MAC-project face velocities, but no divu
      macProject(AmrOldTime,advData,{});

   } else {

      // SDC iterations
      for (int sdc_iter = 1; sdc_iter <= m_nSDCmax; ++sdc_iter ) {
         oneSDC(sdc_iter,advData,diffData);
      }

      // Post SDC
      averageDownDensity(AmrNewTime); // Gather the following if needed TODO
      averageDownSpecies(AmrNewTime);
      averageDownEnthalpy(AmrNewTime);
      averageDownTemp(AmrNewTime);
      fillPatchState(AmrNewTime);
      if (m_has_divu) {
         int is_initialization = 0;             // Not here
         int computeDiffusionTerm = 1;          // Yes, re-evaluate the diffusion term after the last chemistry solve
         int do_avgDown = 1;                    // Always
         calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
      }
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Velocity advance
   Real VelAdvStart = 0.0;
   if (m_verbose > 1) {
      VelAdvStart = ParallelDescriptor::second();
   }
   // Re-evaluate viscosity only if scalar updated
   if (!m_incompressible) calcViscosity(AmrNewTime);

   // Compute t^{n+1/2} velocity advection term
   computeVelocityAdvTerm(advData);

   // Compute provisional new velocity for diffusion solve RHS
   // U^{np1**} = U^{n} - dt*AofS^{n+1/2} - dt/rho^{n+1/2} \nabla \pi^{n-1/2} + dt/rho^{n+1/2} * F^{n+1/2}
   updateVelocity(is_initIter,advData);

   // Semi-implicit CN diffusion solve to get U^{np1*}
   diffuseVelocity();

   // Nodal projection to get constrained U^{np1} and new pressure \pi^{n+1/2}
   const TimeStamp rhoTime = AmrHalfTime;
   velocityProjection(is_initIter,rhoTime,m_dt);
   if (m_verbose > 1) {
      Real VelAdvEnd = ParallelDescriptor::second() - VelAdvStart;
      ParallelDescriptor::ReduceRealMax(VelAdvEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - Advance()::VelocityAdvance " << VelAdvEnd << "\n";
   }
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

   if (m_verbose > 0) {
      amrex::Print() << "   SDC iter [" << sdcIter << "] \n";
   }

   //----------------------------------------------------------------
   // Update t^{n+1,k} transport/Dnp1/divU
   //----------------------------------------------------------------
   // At the first SDC, we already copied old -> new
   if (sdcIter > 1) {

      // fillpatch the new state
      averageDownDensity(AmrNewTime); // Gather the following if needed TODO
      averageDownSpecies(AmrNewTime);
      averageDownEnthalpy(AmrNewTime);
      averageDownTemp(AmrNewTime);
      fillPatchState(AmrNewTime);
      
      calcDiffusivity(AmrNewTime);
      computeDifferentialDiffusionTerms(AmrNewTime,diffData);
      if (m_has_divu) {
         int is_initialization = 0;                   // Not here 
         int computeDiffusionTerm = 0;                // Nope, we just did that
         int do_avgDown = 1;                          // Always
         calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
      }
#ifdef PLM_USE_EFIELD
      ionDriftVelocity(advData);
#endif
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Get u MAC
   //----------------------------------------------------------------
   Real MACStart = 0.0;
   if (m_verbose > 1) {
      MACStart = ParallelDescriptor::second();
   }
   // Predict face velocity with Godunov
   predictVelocity(advData,diffData);

   // Create S^{n+1/2} by fillpatching t^{n} and t^{np1,k}
   createMACRHS(advData);

   // Re-evaluate thermo. pressure and add chi_increment
   addChiIncrement(sdcIter, AmrNewTime,advData);

   // MAC projection
   macProject(AmrOldTime,advData,GetVecOfPtrs(advData->mac_divu));
   if (m_verbose > 1) {
      Real MACEnd = ParallelDescriptor::second() - MACStart;
      ParallelDescriptor::ReduceRealMax(MACEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - oneSDC()::MACProjection() " << MACEnd << "\n";
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar advections
   //----------------------------------------------------------------
   Real ScalAdvStart = 0.0;
   if (m_verbose > 1) {
      ScalAdvStart = ParallelDescriptor::second();
   }
   // Get scalar advection SDC forcing
   getScalarAdvForce(advData,diffData);

   // Get AofS: (\nabla \cdot (\rho Y Umac))^{n+1/2,k}
   // and for density = \sum_k AofS_k
   computeScalarAdvTerms(advData);

   // Compute \rho^{np1,k+1} and fillpatch new density
   updateDensity(advData);
   fillPatchDensity(AmrNewTime);
   if (m_verbose > 1) {
      Real ScalAdvEnd = ParallelDescriptor::second() - ScalAdvStart;
      ParallelDescriptor::ReduceRealMax(ScalAdvEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - oneSDC()::ScalarAdvection() " << ScalAdvEnd << "\n";
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar diffusion
   //----------------------------------------------------------------
   Real ScalDiffStart = 0.0;
   if (m_verbose > 1) {
      ScalDiffStart = ParallelDescriptor::second();
   }
   // Get scalar diffusion SDC RHS (stored in Forcing)
   getScalarDiffForce(advData,diffData);

   // Diffuse scalars
   differentialDiffusionUpdate(advData,diffData);
   if (m_verbose > 1) {
      Real ScalDiffEnd = ParallelDescriptor::second() - ScalDiffStart;
      ParallelDescriptor::ReduceRealMax(ScalDiffEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - oneSDC()::ScalarDiffusion() " << ScalDiffEnd << "\n";
   }
   //----------------------------------------------------------------

#ifdef PLM_USE_EFIELD
   //----------------------------------------------------------------
   // Solve for implicit non-linear nE/PhiV system
   //----------------------------------------------------------------
   implicitNonLinearSolve(sdcIter, m_dt, diffData, advData);
#endif

   //----------------------------------------------------------------
   // Reaction
   //----------------------------------------------------------------
   Real ScalReacStart = 0.0;
   if (m_verbose > 1) {
      ScalReacStart = ParallelDescriptor::second();
   }
   // Get external forcing for chemistry
   getScalarReactForce(advData);

   // Integrate chemistry
   advanceChemistry(advData);
   if (m_verbose > 1) {
      Real ScalReacEnd = ParallelDescriptor::second() - ScalReacStart;
      ParallelDescriptor::ReduceRealMax(ScalReacEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - oneSDC()::ScalarReaction() " << ScalReacEnd << "\n";
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Wrap it up
   //----------------------------------------------------------------
   // Re-evaluate derived state entries
   setTemperature(AmrNewTime);
   floorSpecies(AmrNewTime);
   setThermoPress(AmrNewTime);
}
