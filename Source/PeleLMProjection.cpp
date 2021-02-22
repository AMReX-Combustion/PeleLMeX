#include <PeleLM.H>

using namespace amrex;

void PeleLM::initialProjection()
{
   BL_PROFILE_VAR("PeleLM::initialProjection()", initialProjection);

   Real time = 0.0;

   if (m_verbose) {
      amrex::Print() << "Initial projection: \n";
   }

   Real dummy_dt = 1.0;
   int increment_gp = 0;
   int nGhost = 0;
   
   // Get sigma : density if not incompressible
   Vector<std::unique_ptr<MultiFab>> sigma(finest_level+1);
   if (! m_incompressible ) {
      for (int lev = 0; lev <= finest_level; ++lev ) {
         
         sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

         auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldata_p->density,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& rho_arr = ldata_p->density.const_array(mfi);
            auto const& sig_arr = sigma[lev]->array(mfi);
            amrex::ParallelFor(bx, [rho_arr,sig_arr,dummy_dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               sig_arr(i,j,k) = dummy_dt / rho_arr(i,j,k);
            });
         }
         
      }
   }

   // Get velocity
   Vector<MultiFab*> vel;
   for (int lev = 0; lev <= finest_level; ++lev) {
      vel.push_back(&(m_leveldata_new[lev]->velocity));
      vel[lev]->setBndry(0.0);
      setPhysBoundaryVel(*vel[lev],lev,AmrNewTime);
   }

   // Get RHS cc
   Vector<MultiFab*> rhs_cc;
   if (!m_incompressible && m_has_divu ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         rhs_cc.push_back(&(m_leveldata_new[lev]->divu));
      }
   }

   doNodalProject(vel, amrex::GetVecOfPtrs(sigma), rhs_cc, {}, increment_gp, dummy_dt);

   // Set back press and gpress to zero
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      ldata_p->press.setVal(0.0);
      ldata_p->gp.setVal(0.0);
   }

   if (m_verbose) {
      amrex::Print() << "After initial projection:\n";
   }
  
}

void PeleLM::doNodalProject(Vector<MultiFab*> &a_vel,
                            const Vector<MultiFab*> &a_sigma,
                            const Vector<MultiFab*> &rhs_cc,
                            const Vector<const MultiFab*> &rhs_nd,
                            int increment_gp, 
                            Real scaling_factor) {

   // Asserts
   
   LPInfo info;
   info.setMaxCoarseningLevel(m_nodal_mg_max_coarsening_level);

   int has_rhs = 0;
   if (!rhs_cc.empty()) has_rhs = 1;
   
   // BCs
   std::array<LinOpBCType,AMREX_SPACEDIM> lobc;
   std::array<LinOpBCType,AMREX_SPACEDIM> hibc;
   for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      if (Geom(0).isPeriodic(idim)) {
         lobc[idim] = hibc[idim] = LinOpBCType::Periodic;
      } else {
         if (m_phys_bc.lo(idim) == Outflow) {
            lobc[idim] = LinOpBCType::Dirichlet;
         } else if (m_phys_bc.lo(idim) == Inflow) {
            lobc[idim] = LinOpBCType::inflow;
         } else {
            lobc[idim] = LinOpBCType::Neumann;
         }
         if (m_phys_bc.hi(idim) == Outflow) {
            hibc[idim] = LinOpBCType::Dirichlet;
         } else if (m_phys_bc.hi(idim) == Inflow) {
            hibc[idim] = LinOpBCType::inflow;
         } else {
            hibc[idim] = LinOpBCType::Neumann;
         }
      }
   }

   // Setup NodalProjector
   std::unique_ptr<NodalProjector> nodal_projector;

   if ( m_incompressible ) {
      Real constant_sigma = scaling_factor / m_rho;
      nodal_projector.reset(new NodalProjector(a_vel, constant_sigma, Geom(0,finest_level), info));
   } else {
      if ( has_rhs ) {
         nodal_projector.reset(new NodalProjector(a_vel, GetVecOfConstPtrs(a_sigma), Geom(0,finest_level),
                                                  info, rhs_cc, rhs_nd)); 
      } else {
         nodal_projector.reset(new NodalProjector(a_vel, GetVecOfConstPtrs(a_sigma), Geom(0,finest_level), info)); 
      }
   }
   
   nodal_projector->setDomainBC(lobc, hibc);

#if (AMREX_SPACEDIM == 2)
   if (m_rz_correction) {
      nodal_projector->getLinOp().setRZCorrection(Geom(0).IsRZ());
   }   
#endif

   // Solve
   nodal_projector->project(m_nodal_mg_rtol, m_nodal_mg_atol);

   auto phi = nodal_projector->getPhi();
   auto gphi = nodal_projector->getGradPhi();

   for(int lev = 0; lev <= finest_level; lev++) {

      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef _OPENMP   
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->gp,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
         Box const& tbx = mfi.tilebox();
         Box const& nbx = mfi.nodaltilebox();
         auto const& p_lev_arr   = ldata_p->press.array(mfi);
         auto const& gp_lev_arr  = ldata_p->gp.array(mfi);
         auto const& p_proj_arr  = phi[lev]->const_array(mfi);
         auto const& gp_proj_arr = gphi[lev]->const_array(mfi);
         if (increment_gp) {
            amrex::ParallelFor(tbx, AMREX_SPACEDIM, [gp_lev_arr,gp_proj_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               gp_lev_arr(i,j,k,n) += gp_proj_arr(i,j,k,n);
            });
         } else {
            amrex::ParallelFor(tbx, AMREX_SPACEDIM, [gp_lev_arr,gp_proj_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               gp_lev_arr(i,j,k,n) = gp_proj_arr(i,j,k,n);
            });
         }
         amrex::ParallelFor(nbx, [p_lev_arr,p_proj_arr]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            p_lev_arr(i,j,k) = p_proj_arr(i,j,k);
         });
      }
   }
}
