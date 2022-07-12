#include <PeleLM.H>
#include <PeleLMUtils.H>
#include <AMReX_MemProfiler.H>

using namespace amrex;

void PeleLM::Advance(int is_initIter) {
   BL_PROFILE("PeleLM::Advance()");

#ifdef AMREX_MEM_PROFILING
   // Memory profiler if compiled
   MemProfiler::report("STEP ["+std::to_string(m_nstep)+"]");
#endif

   // Start timing current time step
   Real strt_time = ParallelDescriptor::second();

   //----------------------------------------------------------------
   BL_PROFILE_VAR("PeleLM::advance::setup", PLM_SETUP);
   //----------------------------------------------------------------

   // Deal with ambient pressure
   if (m_closed_chamber) {
      m_pNew = m_pOld;
   }

   // Put together new typical values
   if (!is_initIter && m_nstep > 0 &&
       m_nstep%m_resetTypValInt == 0) {
       setTypicalValues(AmrNewTime);
   }

   //----------------------------------------------------------------
   // Copy old <- new state
   copyStateNewToOld(1);
   copyPressNewToOld();
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // TIME
   // Compute time-step size
   m_dt = computeDt(is_initIter,AmrOldTime);
#ifdef PELELM_USE_SPRAY
   if (!is_initIter && do_spray_particles) {
     // Create the state MF used for spray interpolation
     setSprayState(m_dt);
   }
#endif

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

   checkMemory("Adv. start");

   //----------------------------------------------------------------
   // Data for the advance, only live for the duration of the advance
   std::unique_ptr<AdvanceDiffData> diffData;
   diffData.reset(new AdvanceDiffData(finest_level, grids, dmap, m_factory, m_nGrowAdv, m_use_wbar));
   std::unique_ptr<AdvanceAdvData> advData;
   advData.reset(new AdvanceAdvData(finest_level, grids, dmap, m_factory, m_incompressible,
                                    m_nGrowAdv, m_nGrowMAC));

   for (int lev = 0; lev <= finest_level; lev++) {
     m_extSource[lev]->setVal(0.);
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Advance setup

   // initiliaze temporals
   initTemporals();

   // Compute velocity flux on boundary faces if doing closed chamber
   if (m_closed_chamber) {
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_domainUmacFlux[2*idim] = 0.0;
         m_domainUmacFlux[2*idim+1] = 0.0;
      }
   }

   // fillpatch the t^{n} data
   averageDownState(AmrOldTime);
   fillPatchState(AmrOldTime);

   // compute t^{n} data
   calcViscosity(AmrOldTime);
   if (! m_incompressible ) {
      calcDiffusivity(AmrOldTime);
#ifdef PELE_USE_EFIELD
      poissonSolveEF(AmrOldTime);
#endif
   }
   // TODO : check dt

   //----------------------------------------------------------------
   BL_PROFILE_VAR_STOP(PLM_SETUP);
   //----------------------------------------------------------------

#ifdef PELELM_USE_SPRAY
   if (!is_initIter) {
     sprayMKD(m_cur_time, m_dt);
   }
#endif
#ifdef PELELM_USE_SOOT
   if (do_soot_solve) {
     computeSootSource(AmrOldTime, m_dt);
   }
#endif

   if (! m_incompressible ) {
      floorSpecies(AmrOldTime);

      //----------------------------------------------------------------
      BL_PROFILE_VAR("PeleLM::advance::mac", PLM_MAC);
      setThermoPress(AmrOldTime);
      BL_PROFILE_VAR_STOP(PLM_MAC);
      //----------------------------------------------------------------
      BL_PROFILE_VAR("PeleLM::advance::diffusion", PLM_DIFF);
      computeDifferentialDiffusionTerms(AmrOldTime,diffData);
      BL_PROFILE_VAR_STOP(PLM_DIFF);
      //----------------------------------------------------------------
   }

   // Initialize terms t^{np1,k} from t^{n}
   //----------------------------------------------------------------
   BL_PROFILE_VAR_START(PLM_SETUP);
   //----------------------------------------------------------------
   copyTransportOldToNew();
   if (! m_incompressible ) {
      copyDiffusionOldToNew(diffData);
#ifdef PELE_USE_EFIELD
      ionDriftVelocity(advData);
#endif
   }
   BL_PROFILE_VAR_STOP(PLM_SETUP);

   // TODO : handle reaction ghost cells
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar advance
   if ( m_incompressible ) {
      Real MACStart = 0.0;
      if (m_verbose > 1) {
         MACStart = ParallelDescriptor::second();
      }

      // Still need to get face velocities ...
      predictVelocity(advData);

      // ... and MAC-project face velocities, but no divu
      macProject(AmrOldTime,advData,{});

      if (m_verbose > 1) {
         Real MACEnd = ParallelDescriptor::second() - MACStart;
         ParallelDescriptor::ReduceRealMax(MACEnd, ParallelDescriptor::IOProcessorNumber());
         amrex::Print() << "   - Advance()::MACProjection()  --> Time: " << MACEnd << "\n";
      }
      checkMemory("MAC-Proj");

   } else {

      // SDC iterations
      for (int sdc_iter = 1; sdc_iter <= m_nSDCmax; ++sdc_iter ) {
         oneSDC(sdc_iter,advData,diffData);
      }

      // Post SDC
      averageDownScalars(AmrNewTime);
      fillPatchState(AmrNewTime);

#ifdef PELELM_USE_SOOT
      if (do_soot_solve) {
         clipSootMoments();
      }
#endif
      if (m_has_divu) {
         int is_initialization = 0;             // Not here
         int computeDiffusionTerm = 1;          // Yes, re-evaluate the diffusion term after the last chemistry solve
         int do_avgDown = 1;                    // Always
         calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
      }
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   BL_PROFILE_VAR("PeleLM::advance::velocity", PLM_VEL);
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
   updateVelocity(advData);

   // Semi-implicit CN diffusion solve to get U^{np1*}
   diffuseVelocity();

   // Nodal projection to get constrained U^{np1} and new pressure \pi^{n+1/2}
   const TimeStamp rhoTime = AmrHalfTime;
   velocityProjection(is_initIter,rhoTime,m_dt);
   if (m_verbose > 1) {
      Real VelAdvEnd = ParallelDescriptor::second() - VelAdvStart;
      ParallelDescriptor::ReduceRealMax(VelAdvEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - Advance()::VelocityAdvance  --> Time: " << VelAdvEnd << "\n";
   }
   checkMemory("Nodal-Proj");
   BL_PROFILE_VAR_STOP(PLM_VEL);
   //----------------------------------------------------------------

   // Deal with ambient pressure
   if (m_closed_chamber && !is_initIter) {
      m_pOld = m_pNew;
   }

   //----------------------------------------------------------------
   // Wrapup advance
   // Timing current time step
   if (m_verbose > 0)
   {
      Real run_time = ParallelDescriptor::second() - strt_time;
      ParallelDescriptor::ReduceRealMax(run_time, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << " >> PeleLM::Advance() --> Time: " << run_time << "\n";
   }
}

void PeleLM::oneSDC(int sdcIter,
                    std::unique_ptr<AdvanceAdvData> &advData,
                    std::unique_ptr<AdvanceDiffData> &diffData)
{
   BL_PROFILE("PeleLM::oneSDC()");
   m_sdcIter = sdcIter;

   if (m_verbose > 0) {
      amrex::Print() << "   SDC iter [" << sdcIter << "] \n";
   }

   //----------------------------------------------------------------
   // Update t^{n+1,k} transport/Dnp1/divU
   //----------------------------------------------------------------
   // At the first SDC, we already copied old -> new
   if (sdcIter > 1) {

      Real UpdateStart = 0.0;
      if (m_verbose > 1) {
         UpdateStart = ParallelDescriptor::second();
      }
      // fillpatch the new state
      averageDownScalars(AmrNewTime);
      fillPatchState(AmrNewTime);

      calcDiffusivity(AmrNewTime);
      computeDifferentialDiffusionTerms(AmrNewTime,diffData);
      if (m_has_divu) {
         int is_initialization = 0;                   // Not here
         int computeDiffusionTerm = 0;                // Nope, we just did that
         int do_avgDown = 1;                          // Always
         calcDivU(is_initialization,computeDiffusionTerm,do_avgDown,AmrNewTime,diffData);
      }
#ifdef PELE_USE_EFIELD
      ionDriftVelocity(advData);
#endif
      if (m_verbose > 1) {
         Real UpdateEnd = ParallelDescriptor::second() - UpdateStart;
         ParallelDescriptor::ReduceRealMax(UpdateEnd, ParallelDescriptor::IOProcessorNumber());
         amrex::Print() << "   - oneSDC()::Update t^{n+1,k}  --> Time: " << UpdateEnd << "\n";
      }
   }
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Get u MAC
   //----------------------------------------------------------------
   BL_PROFILE_VAR("PeleLM::advance::mac", PLM_MAC);
   Real MACStart = 0.0;
   if (m_verbose > 1) {
      MACStart = ParallelDescriptor::second();
   }
   // Predict face velocity with Godunov
   predictVelocity(advData);

   // Create S^{n+1/2} by fillpatching t^{n} and t^{np1,k}
   createMACRHS(advData);

   // Re-evaluate thermo. pressure and add chi_increment
   addChiIncrement(sdcIter, AmrNewTime, advData);

   // MAC projection
   macProject(AmrOldTime,advData,GetVecOfPtrs(advData->mac_divu));
   if (m_verbose > 1) {
      Real MACEnd = ParallelDescriptor::second() - MACStart;
      ParallelDescriptor::ReduceRealMax(MACEnd, ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   - oneSDC()::MACProjection()   --> Time: " << MACEnd << "\n";
   }
   checkMemory("MAC-Proj");
   BL_PROFILE_VAR_STOP(PLM_MAC);
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar advections
   //----------------------------------------------------------------
   BL_PROFILE_VAR("PeleLM::advance::scalars_adv", PLM_SADV);
   Real ScalAdvStart = 0.0;
   if (m_verbose > 1) {
      ScalAdvStart = ParallelDescriptor::second();
   }
#ifdef PELELM_USE_SOOT
   // Compute and update passive advective terms
   computePassiveAdvTerms(advData, FIRSTSOOT, NUMSOOTVAR);
#endif
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
      amrex::Print() << "   - oneSDC()::ScalarAdvection() --> Time: " << ScalAdvEnd << "\n";
   }
   checkMemory("ScalAdv");
   BL_PROFILE_VAR_STOP(PLM_SADV);
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Scalar diffusion
   //----------------------------------------------------------------
   BL_PROFILE_VAR("PeleLM::advance::diffusion", PLM_DIFF);
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
      amrex::Print() << "   - oneSDC()::ScalarDiffusion() --> Time: " << ScalDiffEnd << "\n";
   }
   checkMemory("ScalDiff");
   BL_PROFILE_VAR_STOP(PLM_DIFF);
   //----------------------------------------------------------------

#ifdef PELE_USE_EFIELD
   //----------------------------------------------------------------
   // Solve for implicit non-linear nE/PhiV system
   //----------------------------------------------------------------
   implicitNonLinearSolve(sdcIter, m_dt, diffData, advData);
#endif

   //----------------------------------------------------------------
   // Reaction
   //----------------------------------------------------------------
   BL_PROFILE_VAR("PeleLM::advance::reactions", PLM_REAC);
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
      amrex::Print() << "   - oneSDC()::ScalarReaction()  --> Time: " << ScalReacEnd << "\n";
   }
   checkMemory("ScalReact");
   BL_PROFILE_VAR_STOP(PLM_REAC);
   //----------------------------------------------------------------

   //----------------------------------------------------------------
   // Wrap it up
   //----------------------------------------------------------------
   // Re-evaluate derived state entries
   setTemperature(AmrNewTime);
   floorSpecies(AmrNewTime);
   setThermoPress(AmrNewTime);
}
