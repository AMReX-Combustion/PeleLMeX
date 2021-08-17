#include <PeleLM.H>

using namespace amrex;

void PeleLM::initTemporals()
{
   if (!m_do_temporals
       && !(m_nstep % m_temp_int == 0)) return;

   // Reset mass fluxes integrals on domain boundaries
   if (m_do_massBalance) {
      m_massOld = MFSum(GetVecOfConstPtrs(getDensityVect(AmrOldTime)),0);
      for (int idim = 0; idim <  AMREX_SPACEDIM; idim++) {
         m_domainMassFlux[2*idim] = 0.0;
         m_domainMassFlux[2*idim+1] = 0.0;
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

      Real sumLo = 0.0;
      Real sumHi = 0.0;

      sumLo = amrex::ReduceSum(*a_fluxes[idim], 0, [=] 
         AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& flux) -> Real
         {
            Real t = 0.0;
            AMREX_LOOP_3D(bx, i, j, k,
            {
               int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
               if ( idx == faceDomain.smallEnd(idim) ) {
                  for (int n = 0; n < NUM_SPECIES; n++) {
                     t += flux(i,j,k,n) * area[idim];
                  }
               }
            });
            return t;
         });

      sumHi = amrex::ReduceSum(*a_fluxes[idim], 0, [=] 
         AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& flux) -> Real
         {
            Real t = 0.0;
            AMREX_LOOP_3D(bx, i, j, k,
            {
               int idx = (idim == 0) ? i : ( (idim == 1) ? j : k );
               if ( idx == faceDomain.bigEnd(idim) ) {
                  for (int n = 0; n < NUM_SPECIES; n++) {
                     t += flux(i,j,k,n) * area[idim];
                  }
               }
            });
            return t;
         });

      ParallelAllReduce::Sum(sumLo, ParallelContext::CommunicatorSub());
      ParallelAllReduce::Sum(sumHi, ParallelContext::CommunicatorSub());
      m_domainMassFlux[2*idim] += sumLo;
      m_domainMassFlux[2*idim+1] -= sumHi;   // Outflow, negate flux
   }
}

void PeleLM::writeTemporals()
{
   //----------------------------------------------------------------
   // Mass balance
   if (m_do_massBalance) {
      massBalance();
   }

   //----------------------------------------------------------------
   // State
   // Get kinetic energy
   Vector<std::unique_ptr<MultiFab>> kinEnergy(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      kinEnergy[lev] = derive("kinetic_energy", m_cur_time, lev, 0);
   }
   Real kinetic_energy = MFSum(GetVecOfConstPtrs(kinEnergy),0);

   // Get min/max/mean for non-species state components

   tmpStateFile << m_nstep << " " << m_cur_time << " " << kinetic_energy << " \n";
   tmpStateFile.flush();
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
   }
}
