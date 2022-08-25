#include <PeleLM.H>

using namespace amrex;

void PeleLM::initTemporals(const PeleLM::TimeStamp &a_time)
{
   if (!m_do_temporals
       && !(m_nstep % m_temp_int == 0)) return;

   // Reset mass fluxes integrals on domain boundaries
   if (m_do_massBalance && !m_incompressible) {
      m_massOld = MFSum(GetVecOfConstPtrs(getDensityVect(a_time)),0);
      for (int idim = 0; idim <  AMREX_SPACEDIM; idim++) {
         m_domainMassFlux[2*idim] = 0.0;
         m_domainMassFlux[2*idim+1] = 0.0;
      }
   }
   if (m_do_energyBalance && !m_incompressible) {
      m_RhoHOld = MFSum(GetVecOfConstPtrs(getRhoHVect(a_time)),0);
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
         m_domainRhoHFlux[2*idim] = 0.0;
         m_domainRhoHFlux[2*idim+1] = 0.0;
      }
   }

   if (m_do_speciesBalance && !m_incompressible) {
      for (int n = 0; n < NUM_SPECIES; n++){
         m_RhoYOld[n] = MFSum(GetVecOfConstPtrs(getSpeciesVect(a_time)),n);
         for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            m_domainRhoYFlux[2*n*AMREX_SPACEDIM+2*idim] = 0.0;
            m_domainRhoYFlux[1+2*n*AMREX_SPACEDIM+2*idim] = 0.0;
         }
      }
   }

}

void PeleLM::massBalance()
{
   // Compute the mass balance on the computational domain
   m_massNew = MFSum(GetVecOfConstPtrs(getDensityVect(AmrNewTime)),0);
   Real dmdt = (m_massNew - m_massOld) / m_dt;
   Real massFluxBalance = AMREX_D_TERM(  m_domainMassFlux[0] + m_domainMassFlux[1],
                                       + m_domainMassFlux[2] + m_domainMassFlux[3],
                                       + m_domainMassFlux[4] + m_domainMassFlux[5]);

   tmpMassFile << m_nstep << " " << m_cur_time                          // Time info
               << " " << m_massNew                                      // mass
               << " " << dmdt                                           // mass temporal derivative
               << " " << massFluxBalance                                // domain boundaries mass fluxes
               << " " << std::abs(dmdt - massFluxBalance) << " \n";     // balance
   tmpMassFile.flush();
}

void PeleLM::speciesBalance()
{
   // Compute the species rhoY balance on the computational domain
   Array<Real,NUM_SPECIES> dmYdt;
   Array<Real,NUM_SPECIES> massYFluxBalance;
   Array<Real,NUM_SPECIES> rhoYdots;
   for (int n = 0; n < NUM_SPECIES; n++){
       m_RhoYNew[n] = MFSum(GetVecOfConstPtrs(getSpeciesVect(AmrNewTime)),n);
       rhoYdots[n] =  MFSum(GetVecOfConstPtrs(getIRVect()),n);
       dmYdt[n] = (m_RhoYNew[n] - m_RhoYOld[n]) / m_dt;
       massYFluxBalance[n] = AMREX_D_TERM(  m_domainRhoYFlux[2*n*AMREX_SPACEDIM]   + m_domainRhoYFlux[2*n*AMREX_SPACEDIM+1],
                                          + m_domainRhoYFlux[2*n*AMREX_SPACEDIM+2] + m_domainRhoYFlux[2*n*AMREX_SPACEDIM+3],
                                          + m_domainRhoYFlux[2*n*AMREX_SPACEDIM+4] + m_domainRhoYFlux[2*n*AMREX_SPACEDIM+5]);
   }

   tmpSpecFile << m_nstep << " " << m_cur_time;                         // Time info
   for (int n = 0; n < NUM_SPECIES; n++){
       tmpSpecFile << " " << m_RhoYNew[n]                               // mass of Y
                   << " " << dmYdt[n]                                   // mass temporal derivative
                   << " " << massYFluxBalance[n]                        // domain boundaries mass fluxes
                   << " " << rhoYdots[n]                                // integrated consumption rate
                   << " " << std::abs(dmYdt[n] - massYFluxBalance[n]
                                      - rhoYdots[n]);                   // balance
   }
   tmpSpecFile << "\n";
   tmpSpecFile.flush();
}

void PeleLM::addMassFluxes(const Array<const MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                           const Geometry& a_geom)
{

   // Do when m_nstep is -1 since m_nstep is increased by one before
   // the writeTemporals
   if ( !(m_nstep % m_temp_int == m_temp_int-1) ) return;

   // Get the face areas
   const Real*  dx = a_geom.CellSize();
   Array<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 1 )
   area[0] = 1.0;
#elif ( AMREX_SPACEDIM == 2 )
   area[0] = dx[1];
   area[1] = dx[0];
#else
   area[0] = dx[1]*dx[2];
   area[1] = dx[0]*dx[2];
   area[2] = dx[0]*dx[1];
#endif

   for (int idim = 0; idim <  AMREX_SPACEDIM; idim++) {
      auto faceDomain = amrex::convert(a_geom.Domain(),IntVect::TheDimensionVector(idim));

      auto const& fma = a_fluxes[idim]->const_arrays();

      Real sumLo = 0.0;
      Real sumHi = 0.0;

#if ( AMREX_SPACEDIM == 2 )
      if (geom[0].IsRZ()) {
         MultiFab mf_a;
         geom[0].GetFaceArea(mf_a, grids[0], dmap[0], idim, 0); 
         auto const& ama = mf_a.const_arrays();
         auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                   *a_fluxes[idim],IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                                   {
                                       Array4<const Real> const& flux = fma[box_no];
                                       Array4<const Real> const& area_ar = ama[box_no];

                                       int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                       // low
                                       Real low = 0.0;
                                       if ( idx == faceDomain.smallEnd(idim) ) {
                                          for (int n = 0; n < NUM_SPECIES; n++) {
                                             low += flux(i,j,k,n) * area_ar(i,j,k);
                                          }
                                       }
                                       // high
                                       Real high = 0.0;
                                       if ( idx == faceDomain.bigEnd(idim) ) {
                                          for (int n = 0; n < NUM_SPECIES; n++) {
                                             high += flux(i,j,k,n) * area_ar(i,j,k);
                                          }
                                       }
                                       return{low,high};
                                   });
         sumLo = amrex::get<0>(r);
         sumHi = amrex::get<1>(r);
      } else
#endif
      {
         auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                   *a_fluxes[idim],IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                                   {
                                       Array4<const Real> const& flux = fma[box_no];

                                       int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                       // low
                                       Real low = 0.0;
                                       if ( idx == faceDomain.smallEnd(idim) ) {
                                          for (int n = 0; n < NUM_SPECIES; n++) {
                                             low += flux(i,j,k,n) * area[idim];
                                          }
                                       }
                                       // high
                                       Real high = 0.0;
                                       if ( idx == faceDomain.bigEnd(idim) ) {
                                          for (int n = 0; n < NUM_SPECIES; n++) {
                                             high += flux(i,j,k,n) * area[idim];
                                          }
                                       }
                                       return{low,high};
                                   });
         sumLo = amrex::get<0>(r);
         sumHi = amrex::get<1>(r);
      }
      ParallelAllReduce::Sum<Real>({sumLo,sumHi}, ParallelContext::CommunicatorSub());
      m_domainMassFlux[2*idim] += sumLo;
      m_domainMassFlux[2*idim+1] -= sumHi;   // Outflow, negate flux
   }
}

void PeleLM::addUmacFluxes(std::unique_ptr<AdvanceAdvData> &advData, const Geometry& a_geom)
{
   // Get the face areas
   const Real*  dx = a_geom.CellSize();
   Array<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 1 )
   area[0] = 1.0;
#elif ( AMREX_SPACEDIM == 2 )
   area[0] = dx[1];
   area[1] = dx[0];
#else
   area[0] = dx[1]*dx[2];
   area[1] = dx[0]*dx[2];
   area[2] = dx[0]*dx[1];
#endif

   // Just use level 0 since we are calling after averaging down
   int lev = 0;

   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      auto faceDomain = amrex::convert(a_geom.Domain(),IntVect::TheDimensionVector(idim));

      auto const& fma = advData->umac[lev][idim].const_arrays();

      Real sumLo = 0.0;
      Real sumHi = 0.0;

#if ( AMREX_SPACEDIM == 2 )
      if (geom[0].IsRZ()) {
         MultiFab mf_a;
         geom[0].GetFaceArea(mf_a, grids[0], dmap[0], idim, 0); 
         auto const& ama = mf_a.const_arrays();
         auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                   advData->umac[lev][idim],IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                                   {
                                       Array4<const Real> const& flux = fma[box_no];
                                       Array4<const Real> const& area_ar = ama[box_no];

                                       int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                       // low
                                       Real low = 0.0;
                                       if ( idx == faceDomain.smallEnd(idim) ) {
                                          low += flux(i,j,k) * area_ar(i,j,k);
                                       }
                                       // high
                                       Real high = 0.0;
                                       if ( idx == faceDomain.bigEnd(idim) ) {
                                          high += flux(i,j,k) * area_ar(i,j,k);
                                       }
                                       return{low,high};
                                   });
         sumLo = amrex::get<0>(r);
         sumHi = amrex::get<1>(r);
      } else
#endif
      {
         auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                   advData->umac[lev][idim],IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                                   {
                                       Array4<const Real> const& flux = fma[box_no];

                                       int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                       // low
                                       Real low = 0.0;
                                       if ( idx == faceDomain.smallEnd(idim) ) {
                                          low += flux(i,j,k) * area[idim];
                                       }
                                       // high
                                       Real high = 0.0;
                                       if ( idx == faceDomain.bigEnd(idim) ) {
                                          high += flux(i,j,k) * area[idim];
                                       }
                                       return{low,high};
                                   });
         sumLo = amrex::get<0>(r);
         sumHi = amrex::get<1>(r);
      }
      ParallelAllReduce::Sum<Real>({sumLo,sumHi}, ParallelContext::CommunicatorSub());
      m_domainUmacFlux[2*idim] += sumLo;
      m_domainUmacFlux[2*idim+1] -= sumHi;   // Outflow, negate flux
   }
}

void PeleLM::rhoHBalance()
{
   // Compute the enthalpy balance on the computational domain (rho*h)
   m_RhoHNew = MFSum(GetVecOfConstPtrs(getRhoHVect(AmrNewTime)),0);
   Real dRhoHdt = (m_RhoHNew - m_RhoHOld) / m_dt;
   Real rhoHFluxBalance = AMREX_D_TERM(  m_domainRhoHFlux[0] + m_domainRhoHFlux[1],
                                       + m_domainRhoHFlux[2] + m_domainRhoHFlux[3],
                                       + m_domainRhoHFlux[4] + m_domainRhoHFlux[5]);

   tmpMassFile << m_nstep << " " << m_cur_time                          // Time info
               << " " << m_RhoHNew                                      // RhoH
               << " " << dRhoHdt                                        // RhoH temporal derivative
               << " " << rhoHFluxBalance                                // domain boundaries RhoH fluxes
               << " " << std::abs(dRhoHdt - rhoHFluxBalance) << " \n";  // balance
   tmpMassFile.flush();
}

void PeleLM::addRhoHFluxes(const Array<const MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                           const Geometry& a_geom)
{

   // Do when m_nstep is -1 since m_nstep is increased by one before
   // the writeTemporals
   if ( !(m_nstep % m_temp_int == m_temp_int-1) ) return;

   // Get the face areas
   const Real*  dx = a_geom.CellSize();
   Array<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 1 )
   area[0] = 1.0;
#elif ( AMREX_SPACEDIM == 2 )
   area[0] = dx[1];
   area[1] = dx[0];
#else
   area[0] = dx[1]*dx[2];
   area[1] = dx[0]*dx[2];
   area[2] = dx[0]*dx[1];
#endif

   for (int idim = 0; idim <  AMREX_SPACEDIM; idim++) {
      auto faceDomain = amrex::convert(a_geom.Domain(),IntVect::TheDimensionVector(idim));

      auto const& fma = a_fluxes[idim]->const_arrays();

      Real sumLo = 0.0;
      Real sumHi = 0.0;

#if ( AMREX_SPACEDIM == 2 )
      if (geom[0].IsRZ()) {
         MultiFab mf_a;
         geom[0].GetFaceArea(mf_a, grids[0], dmap[0], idim, 0); 
         auto const& ama = mf_a.const_arrays();
         auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                   *a_fluxes[idim],IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                                   {
                                       Array4<const Real> const& flux = fma[box_no];
                                       Array4<const Real> const& area_ar = ama[box_no];

                                       int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                       // low
                                       Real low = 0.0;
                                       if ( idx == faceDomain.smallEnd(idim) ) {
                                          low += flux(i,j,k,NUM_SPECIES) * area_ar(i,j,k);
                                       }
                                       // high
                                       Real high = 0.0;
                                       if ( idx == faceDomain.bigEnd(idim) ) {
                                          high += flux(i,j,k,NUM_SPECIES) * area_ar(i,j,k);
                                       }
                                       return{low,high};
                                   });
         sumLo = amrex::get<0>(r);
         sumHi = amrex::get<1>(r);
      } else
#endif
      {
         auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                   *a_fluxes[idim],IntVect(0),
                                   [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                   -> GpuTuple<Real,Real>
                                   {
                                       Array4<const Real> const& flux = fma[box_no];

                                       int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                       // low
                                       Real low = 0.0;
                                       if ( idx == faceDomain.smallEnd(idim) ) {
                                          low += flux(i,j,k,NUM_SPECIES) * area[idim];
                                       }
                                       // high
                                       Real high = 0.0;
                                       if ( idx == faceDomain.bigEnd(idim) ) {
                                          high += flux(i,j,k,NUM_SPECIES) * area[idim];
                                       }
                                       return{low,high};
                                   });
         sumLo = amrex::get<0>(r);
         sumHi = amrex::get<1>(r);
      }
      ParallelAllReduce::Sum<Real>({sumLo,sumHi}, ParallelContext::CommunicatorSub());
      m_domainRhoHFlux[2*idim] += sumLo;
      m_domainRhoHFlux[2*idim+1] -= sumHi;   // Outflow, negate flux
   }
}


void PeleLM::addRhoYFluxes(const Array<const MultiFab*,AMREX_SPACEDIM> &a_fluxes,
                           const Geometry& a_geom, const Real &a_factor)
{

   // Do when m_nstep is -1 since m_nstep is increased by one before
   // the writeTemporals
   if ( !(m_nstep % m_temp_int == m_temp_int-1) ) return;

   // Get the face areas
   const Real*  dx = a_geom.CellSize();
   Array<Real,AMREX_SPACEDIM> area;
#if ( AMREX_SPACEDIM == 1 )
   area[0] = 1.0;
#elif ( AMREX_SPACEDIM == 2 )
   area[0] = dx[1];
   area[1] = dx[0];
#else
   area[0] = dx[1]*dx[2];
   area[1] = dx[0]*dx[2];
   area[2] = dx[0]*dx[1];
#endif

   // Outer loop over species
   for (int n = 0; n < NUM_SPECIES; n++){
      // Inner loop over dimensions
      for (int idim = 0; idim <  AMREX_SPACEDIM; idim++) {
         auto faceDomain = amrex::convert(a_geom.Domain(),IntVect::TheDimensionVector(idim));

         auto const& fma = a_fluxes[idim]->const_arrays();

         Real sumLo = 0.0;
         Real sumHi = 0.0;

#if ( AMREX_SPACEDIM == 2 )
         if (geom[0].IsRZ()) {
            MultiFab mf_a;
            geom[0].GetFaceArea(mf_a, grids[0], dmap[0], idim, 0); 
            auto const& ama = mf_a.const_arrays();
            auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                      *a_fluxes[idim],IntVect(0),
                                      [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                      -> GpuTuple<Real,Real>
                                      {
                                          Array4<const Real> const& flux = fma[box_no];
                                          Array4<const Real> const& area_ar = ama[box_no];

                                          int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                          // low
                                          Real low = 0.0;
                                          if ( idx == faceDomain.smallEnd(idim) ) {
                                             low += flux(i,j,k,n) * area_ar(i,j,k);
                                          }
                                          // high
                                          Real high = 0.0;
                                          if ( idx == faceDomain.bigEnd(idim) ) {
                                             high += flux(i,j,k,n) * area_ar(i,j,k);
                                          }
                                          return{low,high};
                                      });
            sumLo = amrex::get<0>(r);
            sumHi = amrex::get<1>(r);
         } else
#endif
         {
            auto r = amrex::ParReduce(TypeList<ReduceOpSum,ReduceOpSum>{}, TypeList<Real,Real>{},
                                      *a_fluxes[idim],IntVect(0),
                                      [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
                                      -> GpuTuple<Real,Real>
                                      {
                                          Array4<const Real> const& flux = fma[box_no];

                                          int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
                                          // low
                                          Real low = 0.0;
                                          if ( idx == faceDomain.smallEnd(idim) ) {
                                             low += flux(i,j,k,n) * area[idim];
                                          }
                                          // high
                                          Real high = 0.0;
                                          if ( idx == faceDomain.bigEnd(idim) ) {
                                             high += flux(i,j,k,n) * area[idim];
                                          }
                                          return{low,high};
                                      });
            sumLo = amrex::get<0>(r);
            sumHi = amrex::get<1>(r);
         }
         ParallelAllReduce::Sum<Real>({sumLo,sumHi}, ParallelContext::CommunicatorSub());
         m_domainRhoYFlux[2*idim+n*2*AMREX_SPACEDIM] += a_factor * sumLo;
         m_domainRhoYFlux[2*idim+n*2*AMREX_SPACEDIM+1] -= a_factor * sumHi;   // Outflow, negate flux
      }
   }
}

void PeleLM::writeTemporals()
{
   //----------------------------------------------------------------
   // Mass balance
   if (m_do_massBalance && !m_incompressible) {
      massBalance();
   }

   //----------------------------------------------------------------
   // Species balance
   if (m_do_speciesBalance && !m_incompressible) {
      speciesBalance();
   }

   //----------------------------------------------------------------
   // State
   // Get kinetic energy and enstrophy
   Vector<std::unique_ptr<MultiFab>> kinEnergy(finest_level + 1);
   Vector<std::unique_ptr<MultiFab>> enstrophy(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      kinEnergy[lev] = derive("kinetic_energy", m_cur_time, lev, 0);
      enstrophy[lev] = derive("enstrophy", m_cur_time, lev, 0);
   }
   Real kinenergy_int = MFSum(GetVecOfConstPtrs(kinEnergy),0);
   Real enstrophy_int = MFSum(GetVecOfConstPtrs(enstrophy),0);

   // Combustion
   Real fuelConsumptionInt = 0.0;
   Real heatReleaseRateInt = 0.0;
   if (fuelID > 0 && !(m_chem_integrator == "ReactorNull")) {
       fuelConsumptionInt =  MFSum(GetVecOfConstPtrs(getIRVect()),fuelID);
       for (int lev = 0; lev <= finest_level; ++lev) {
          getHeatRelease(lev, kinEnergy[lev].get());  // Re-use kinEnergy container
       }
       heatReleaseRateInt = MFSum(GetVecOfConstPtrs(kinEnergy),0);
   }

   tmpStateFile << m_nstep << " " << m_cur_time << " " << m_dt  // Time
                << " " << kinenergy_int                         // Kinetic energy
                << " " << enstrophy_int                         // Enstrophy
                << " " << m_pNew                                // Thermo. pressure
                << " " << fuelConsumptionInt                    // Integ fuel burning rate
                << " " << heatReleaseRateInt                    // Integ heat release rate
                << " \n";
   tmpStateFile.flush();

   // Get min/max for state components
   auto stateMax  = ( m_incompressible ) ? MLmax(GetVecOfConstPtrs(getStateVect(AmrNewTime)),0,AMREX_SPACEDIM)
                                         : MLmax(GetVecOfConstPtrs(getStateVect(AmrNewTime)),0,NVAR);
   auto stateMin  = ( m_incompressible ) ? MLmin(GetVecOfConstPtrs(getStateVect(AmrNewTime)),0,AMREX_SPACEDIM)
                                         : MLmin(GetVecOfConstPtrs(getStateVect(AmrNewTime)),0,NVAR);

   tmpExtremasFile << m_nstep << " " << m_cur_time;           // Time
   for (int n = 0; n < stateMax.size(); ++n) {                // Min & max of each state variable
       tmpExtremasFile << " " << stateMin[n] << " " << stateMax[n];
   }
   tmpExtremasFile << " \n";
   tmpExtremasFile.flush();


#ifdef PELE_USE_EFIELD
   if (m_do_ionsBalance) {
      ionsBalance();
   }
#endif

}

void PeleLM::openTempFile()
{
   if (!m_do_temporals) return;

   // Create the temporal directory
   UtilCreateDirectory("temporals", 0755);

   if (ParallelDescriptor::IOProcessor()) {
      std::string tempFileName = "temporals/tempState";
      tmpStateFile.open(tempFileName.c_str(),std::ios::out | std::ios::app | std::ios_base::binary);
      tmpStateFile.precision(12);
      if (m_do_massBalance) {
         tempFileName = "temporals/tempMass";
         tmpMassFile.open(tempFileName.c_str(),std::ios::out | std::ios::app | std::ios_base::binary);
         tmpMassFile.precision(12);
      }
      if (m_do_speciesBalance) {
         tempFileName = "temporals/tempSpecies";
         tmpSpecFile.open(tempFileName.c_str(),std::ios::out | std::ios::app | std::ios_base::binary);
         tmpSpecFile.precision(12);
      }
      if (m_do_extremas) {
         tempFileName = "temporals/tempExtremas";
         tmpExtremasFile.open(tempFileName.c_str(),std::ios::out | std::ios::app | std::ios_base::binary);
         tmpExtremasFile.precision(12);
      }
#ifdef PELE_USE_EFIELD
      if (m_do_ionsBalance) {
         tempFileName = "temporals/tempIons";
         tmpIonsFile.open(tempFileName.c_str(),std::ios::out | std::ios::app | std::ios_base::binary);
         tmpIonsFile.precision(12);
      }
#endif
   }
}

void PeleLM::closeTempFile()
{
   if (!m_do_temporals) return;

   if (ParallelDescriptor::IOProcessor()) {
      tmpStateFile.flush();
      tmpStateFile.close();
      if (m_do_massBalance) {
         tmpMassFile.flush();
         tmpMassFile.close();
      }
      if (m_do_speciesBalance) {
         tmpSpecFile.flush();
         tmpSpecFile.close();
      }
      if (m_do_extremas) {
         tmpExtremasFile.flush();
         tmpExtremasFile.close();
      }
#ifdef PELE_USE_EFIELD
      if (m_do_ionsBalance) {
         tmpIonsFile.flush();
         tmpIonsFile.close();
      }
#endif
   }
}
