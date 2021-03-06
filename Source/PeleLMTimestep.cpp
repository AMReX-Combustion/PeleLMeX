#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

Real
PeleLM::computeDt(int is_init,
                  const TimeStamp &a_time)
{
   BL_PROFILE_VAR("PeleLM::computeDt()", computeDt);

   Real estdt = 1.0e200;

   //----------------------------------------------------------------
   // Store prev dt(s)
   m_prev_dt = m_dt;

   //----------------------------------------------------------------
   // Compute dt estimate from level data
   if ( m_fixed_dt > 0.0 ) {
      estdt = m_fixed_dt;
   } else{
      Real dtconv = estConvectiveDt(a_time);
      estdt = std::min(estdt,dtconv);
      if (!m_incompressible && m_has_divu) {
         Real dtdivU = estDivUDt(a_time);
         estdt = std::min(estdt,dtdivU);
#ifdef PLM_USE_EFIELD
         Real dtions = estEFIonsDt(a_time);
         estdt = std::min(estdt, dtions);
#endif
         if ( m_verbose ) {
            Print() << " Est. time step - Conv: " << dtconv << ", divu: " << dtdivU
#ifdef PLM_USE_EFIELD
                    << ", ions: " << dtions
#endif
                    << "\n";
         }
      }
   }

   //----------------------------------------------------------------
   // Limit dt
   if (is_init || m_nstep == 0) {
      estdt *= m_dtshrink;
   } else {
      estdt = std::min(estdt,m_prev_dt*m_dtChangeMax);
      if (m_stop_time >= 0.0) {
         // Ensure ~O(dt) last step by checking a little in advance
         Real timeLeft = (m_stop_time-m_cur_time);
         if ( 2.0 * estdt > timeLeft && timeLeft > estdt ) {
            estdt = 0.5 * timeLeft;
         } else {
            estdt = std::min(estdt,timeLeft);
         }
      }
   }

   return estdt;
}

Real
PeleLM::estConvectiveDt(const TimeStamp &a_time) {

   Real estdt = 1.0e200;
   constexpr Real  small = 1.0e-8;

   for (int lev = 0; lev <= finest_level; ++lev) {

      Real estdt_lev = 1.0e200;

      //----------------------------------------------------------------
      // Get level data ptr
      auto ldata_p = getLevelDataPtr(lev, a_time);
      Real time = getTime(lev, a_time);

      auto const dx = geom[lev].CellSizeArray();

      //----------------------------------------------------------------
      // Get velocity forces
      int nGrow_force = 0;
      MultiFab velForces(grids[lev],dmap[lev],AMREX_SPACEDIM,nGrow_force,MFInfo(),Factory(lev));

      int add_gradP = 1;
      getVelForces(a_time, lev, nullptr, &velForces, add_gradP);

      //----------------------------------------------------------------
      // Get max forces
      Vector<Real> f_max(AMREX_SPACEDIM);
      f_max = velForces.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

      // Get max velocity
      Vector<Real> u_max(AMREX_SPACEDIM);
      u_max = ldata_p->velocity.norm0({AMREX_D_DECL(0,1,2)},0,true,true);

      //----------------------------------------------------------------
      // Est. min time step on lev
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
         if (u_max[idim] > small) {
            estdt_lev = std::min(estdt_lev, dx[idim]/u_max[idim]);
         }
         if (f_max[idim] > small) {
            estdt_lev = std::min(estdt_lev, std::sqrt(2.0*dx[idim]/f_max[idim]));
         }
      }

      //----------------------------------------------------------------
      // Set overall convective dt
      estdt = std::min(estdt,estdt_lev* m_cfl);
   }

   ParallelDescriptor::ReduceRealMin(estdt);

   return estdt;
}

Real
PeleLM::estDivUDt(const TimeStamp &a_time) {

   Real estdt = 1.0e200;

   // Note: only method 1 of PeleLM is available here
   AMREX_ASSERT(m_divu_checkFlag==1);
  
   for (int lev = 0; lev <= finest_level; ++lev) {

      auto ldata_p = getLevelDataPtr(lev, a_time);

      Real divu_dt = amrex::ReduceMin(ldata_p->density, ldata_p->divu, 0,
                                      [dtfac = m_divu_dtFactor, rhoMin = m_divu_rhoMin]
      AMREX_GPU_HOST_DEVICE (Box const& bx, Array4<Real const> const& rho,
                                            Array4<Real const> const& divu ) noexcept -> Real
      {
         using namespace amrex::literals;
         const auto lo = amrex::lbound(bx);
         const auto hi = amrex::ubound(bx);
         amrex::Real dt = 1.e37_rt;
         for       (int k = lo.z; k <= hi.z; ++k) {
            for    (int j = lo.y; j <= hi.y; ++j) {
               for (int i = lo.x; i <= hi.x; ++i) {
                  Real dtcell = est_divu_dt_1(i, j, k, dtfac, rhoMin, rho, divu );
                  dt = amrex::min(dt,dtcell);
               }
            }
         }
         return dt;
      });
      estdt = std::min(divu_dt,estdt);
   }

   ParallelDescriptor::ReduceRealMin(estdt);

   return estdt;
}
