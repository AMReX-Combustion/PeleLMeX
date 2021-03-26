#include <PeleLM.H>

using namespace amrex;

void PeleLM::initialProjection()
{
   BL_PROFILE_VAR("PeleLM::initialProjection()", initialProjection);

   Real time = 0.0;

   if (m_verbose) {
      amrex::Print() << " Initial velocity projection: \n";
   }

   Real dummy_dt = 1.0;
   int incremental = 0;
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

   // Get RHS cc: - divU
   Vector<MultiFab*> rhs_cc;
   if (!m_incompressible && m_has_divu ) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         rhs_cc.push_back(&(m_leveldata_new[lev]->divu));
         rhs_cc[lev]->mult(-1.0,0,1,rhs_cc[lev]->nGrow());
      }
   }

   doNodalProject(vel, amrex::GetVecOfPtrs(sigma), rhs_cc, {}, incremental, dummy_dt);

   // Set back press and gpress to zero and invert divu sign
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      ldata_p->press.setVal(0.0);
      ldata_p->gp.setVal(0.0);
      if (!m_incompressible && m_has_divu ) {
         m_leveldata_new[lev]->divu.mult(-1.0,0,1,rhs_cc[lev]->nGrow());
      }
   }

   // Average down velocity
   for (int lev = finest_level-1; lev >= 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev+1,AmrNewTime);
      auto ldataCrse_p = getLevelDataPtr(lev,AmrNewTime);
#ifdef AMREX_USE_EB
      amrex::EB_average_down(ldataFine_p->velocity, ldataCrse_p->velocity,
                             0, AMREX_SPACEDIM, refRatio(lev));
#else
      amrex::average_down(ldataFine_p->velocity, ldataCrse_p->velocity,
                          0, AMREX_SPACEDIM, refRatio(lev));
#endif
   }

   if (m_verbose) {
      amrex::Print() << " After initial velocity projection:\n";
   }

}

void PeleLM::velocityProjection(int is_initIter,
                                const TimeStamp &a_rhoTime,
                                const Real &a_dt)
{
   BL_PROFILE_VAR("PeleLM::velocityProjection()", velocityProjection);

   int nGhost = 0;
   int incremental = (is_initIter) ? 1 : 0;

   // Get sigma : scaled density inv. if not incompressible
   Vector<std::unique_ptr<MultiFab>> sigma(finest_level+1);
   if (! m_incompressible ) {
      Vector<MultiFab*> rhoHalf(finest_level+1);
      rhoHalf = getDensityVect(a_rhoTime);
      for (int lev = 0; lev <= finest_level; ++lev ) {

         sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(*rhoHalf[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& rho_arr = rhoHalf[lev]->const_array(mfi);
            auto const& sig_arr = sigma[lev]->array(mfi);
            amrex::ParallelFor(bx, [rho_arr,sig_arr,a_dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               sig_arr(i,j,k) = a_dt / rho_arr(i,j,k);
            });
         }
      }
   }

   if (!incremental) {
      Vector<MultiFab*> rhoHalf(finest_level+1);
      if (!m_incompressible) rhoHalf = getDensityVect(a_rhoTime);
      for (int lev = 0; lev <= finest_level; ++lev ) {

         auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
         auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldataNew_p->velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& vel_arr = ldataNew_p->velocity.array(mfi);
            auto const& gp_arr  = ldataOld_p->gp.const_array(mfi);
            auto const& rho_arr = (m_incompressible) ? Array4<Real const>() : rhoHalf[lev]->const_array(mfi);
            amrex::ParallelFor(bx, [vel_arr,gp_arr,rho_arr,a_dt,
                                    incompressible=m_incompressible,rho=m_rho]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               Real soverrho = (incompressible) ? a_dt / rho
                                                : a_dt / rho_arr(i,j,k);
               AMREX_D_TERM(vel_arr(i,j,k,0) += gp_arr(i,j,k,0) * soverrho;,
                            vel_arr(i,j,k,1) += gp_arr(i,j,k,1) * soverrho;,
                            vel_arr(i,j,k,2) += gp_arr(i,j,k,2) * soverrho);
            });
         }
      }
   }

   // If incremental
   // define "vel" to be U^{np1*} - U^{n} rather than U^{np1*}
   if (incremental) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
         auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
         MultiFab::Subtract(ldataNew_p->velocity,
                            ldataOld_p->velocity,
                            0,0,AMREX_SPACEDIM,0);
      }
   }

   // Get velocity
   Vector<MultiFab*> vel;
   for (int lev = 0; lev <= finest_level; ++lev) {
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      vel.push_back(&(ldata_p->velocity));
      vel[lev]->setBndry(0.0);
      if (!incremental) setPhysBoundaryVel(*vel[lev],lev,AmrNewTime);
   }

   // Get RHS cc
   Vector<MultiFab> rhs_cc;
   if (!m_incompressible && m_has_divu ) {
      rhs_cc.resize(finest_level+1);
      for (int lev = 0; lev <= finest_level; ++lev) {
         if (!incremental) {
            auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
            rhs_cc[lev].define(grids[lev],dmap[lev],1,ldata_p->divu.nGrow(), MFInfo(), *m_factory[lev]);
            MultiFab::Copy(rhs_cc[lev],ldata_p->divu,0,0,1,ldata_p->divu.nGrow());
            rhs_cc[lev].mult(-1.0,0,1,ldata_p->divu.nGrow());
         } else {
            auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
            auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
            rhs_cc[lev].define(grids[lev],dmap[lev],1,ldataOld_p->divu.nGrow(), MFInfo(), *m_factory[lev]);
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(rhs_cc[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               const Box& gbx       = mfi.growntilebox();
               const auto& divu_o   = ldataOld_p->divu.const_array(mfi);
               const auto& divu_n   = ldataNew_p->divu.const_array(mfi);
               const auto& rhs      = rhs_cc[lev].array(mfi);
               amrex::ParallelFor(gbx, [divu_o,divu_n,rhs]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  rhs(i,j,k) = - (divu_n(i,j,k) - divu_o(i,j,k));
               });
            }
         }
      }
   }

   doNodalProject(vel, GetVecOfPtrs(sigma), GetVecOfPtrs(rhs_cc), {}, incremental, a_dt);

   // If incremental
   // define back to be U^{np1} by adding U^{n}
   if (incremental) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
         auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
         MultiFab::Add(ldataNew_p->velocity,
                       ldataOld_p->velocity,
                       0,0,AMREX_SPACEDIM,0);
      }
   }

}

void PeleLM::doNodalProject(Vector<MultiFab*> &a_vel,
                            const Vector<MultiFab*> &a_sigma,
                            const Vector<MultiFab*> &rhs_cc,
                            const Vector<const MultiFab*> &rhs_nd,
                            int incremental,
                            Real scaling_factor) {

   // Asserts TODO

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
         if (incremental) {
            amrex::ParallelFor(tbx, AMREX_SPACEDIM, [gp_lev_arr,gp_proj_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               gp_lev_arr(i,j,k,n) += gp_proj_arr(i,j,k,n);
            });
            amrex::ParallelFor(nbx, [p_lev_arr,p_proj_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               p_lev_arr(i,j,k) += p_proj_arr(i,j,k);
            });
         } else {
            amrex::ParallelFor(tbx, AMREX_SPACEDIM, [gp_lev_arr,gp_proj_arr]
            AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
            {
               gp_lev_arr(i,j,k,n) = gp_proj_arr(i,j,k,n);
            });
            amrex::ParallelFor(nbx, [p_lev_arr,p_proj_arr]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               p_lev_arr(i,j,k) = p_proj_arr(i,j,k);
            });
         }
      }
   }


   // Average down grap P
   for (int lev = finest_level-1; lev >= 0; --lev) {
      auto ldataFine_p = getLevelDataPtr(lev+1,AmrNewTime);
      auto ldataCrse_p = getLevelDataPtr(lev,AmrNewTime);
#ifdef AMREX_USE_EB
      amrex::EB_average_down(ldataFine_p->gp, ldataCrse_p->gp,
                             0, AMREX_SPACEDIM, refRatio(lev));
#else
      amrex::average_down(ldataFine_p->gp, ldataCrse_p->gp,
                          0, AMREX_SPACEDIM, refRatio(lev));
#endif
   }

}
