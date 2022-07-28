#include <PeleLM.H>

using namespace amrex;

void PeleLM::initialProjection()
{
   BL_PROFILE_VAR("PeleLM::initialProjection()", initialProjection);

   if (m_verbose) {
      Vector<Real> velMax(AMREX_SPACEDIM);
      velMax = MLNorm0(GetVecOfConstPtrs(getVelocityVect(AmrNewTime)),0,AMREX_SPACEDIM);
      amrex::Print() << " Initial velocity projection: ";
      amrex::Print() << AMREX_D_TERM( "  U: " << velMax[0] <<,
                                      "  V: " << velMax[1] <<,
                                      "  W: " << velMax[2] <<) "\n";
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

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& rho_arr = ldata_p->state.const_array(mfi,DENSITY);
            auto const& sig_arr = sigma[lev]->array(mfi);
            amrex::ParallelFor(bx, [rho_arr,sig_arr,dummy_dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               sig_arr(i,j,k) = dummy_dt / rho_arr(i,j,k);
            });
         }
         scaleProj_RZ(lev,*sigma[lev]);
      }
   }

   // Get velocity
   Vector<std::unique_ptr<MultiFab> > vel;
   for (int lev = 0; lev <= finest_level; ++lev) {
      vel.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state, amrex::make_alias,VELX,AMREX_SPACEDIM));
      vel[lev]->setBndry(0.0);
      setInflowBoundaryVel(*vel[lev],lev,AmrNewTime);
      scaleProj_RZ(lev,*vel[lev]);
   }

   // Get RHS cc: - divU (- \int{divU})
   Real Sbar = 0.0;
   Vector<MultiFab> rhs_cc(finest_level+1);
   if (!m_incompressible && m_has_divu ) {
      // Ensure integral of RHS is zero for closed chamber
      if (m_closed_chamber) {
         Sbar = MFSum(GetVecOfConstPtrs(getDivUVect(AmrNewTime)),0);
         Sbar /= m_uncoveredVol;        // Transform in Mean.
      }
      for (int lev = 0; lev <= finest_level; ++lev) {
         rhs_cc[lev].define(grids[lev],dmap[lev],1,m_leveldata_new[lev]->divu.nGrow());
         MultiFab::Copy(rhs_cc[lev],m_leveldata_new[lev]->divu,0,0,1,m_leveldata_new[lev]->divu.nGrow());
         if (m_closed_chamber) {
            rhs_cc[lev].plus(-Sbar,0,1);
         }
         scaleProj_RZ(lev,rhs_cc[lev]);
         rhs_cc[lev].mult(-1.0,0,1,rhs_cc[lev].nGrow());
      }
   }

   doNodalProject(GetVecOfPtrs(vel), GetVecOfPtrs(sigma), GetVecOfPtrs(rhs_cc), {}, incremental, dummy_dt);

   // Set back press and gpress to zero and restore divu
   // and rescale velocity if 2D-RZ
   for (int lev = 0; lev <= finest_level; lev++) {
      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);
      ldata_p->press.setVal(0.0);
      ldata_p->gp.setVal(0.0);
      if (!m_incompressible && m_has_divu ) {
         m_leveldata_new[lev]->divu.mult(-1.0,0,1,rhs_cc[lev].nGrow());
         // Restore divU integral
         if (m_closed_chamber) {
            m_leveldata_new[lev]->divu.plus(Sbar,0,1);
         }
      }
      unscaleProj_RZ(lev,*vel[lev]);
   }

   // In R-Z, AMReX-Hydro do an average down of r*vel.
   // Now that we have unscaled vel, need to do average down again
   // to have consistent vel across levels
   if (Geom(0).IsRZ()) {
      averageDownVelocity(AmrNewTime);
   }

   if (m_verbose) {
      Vector<Real> velMax(AMREX_SPACEDIM);
      velMax = MLNorm0(GetVecOfConstPtrs(getVelocityVect(AmrNewTime)),0,AMREX_SPACEDIM);
      amrex::Print() << " >> After initial velocity projection: ";
      amrex::Print() << AMREX_D_TERM( "  U: " << velMax[0] <<,
                                      "  V: " << velMax[1] <<,
                                      "  W: " << velMax[2] <<) "\n";
   }
}

void PeleLM::initialPressProjection()
{
   BL_PROFILE_VAR("PeleLM::initialPressProjection()", initialProjection);

   if (m_verbose) {
      amrex::Print() << " Initial pressure projection \n";
   }

   Real dummy_dt = 1.0;
   int incremental = 0;
   int nGhost = 1;

   // Get sigma : density if not incompressible
   Vector<std::unique_ptr<MultiFab>> sigma(finest_level+1);
   if (! m_incompressible ) {
      for (int lev = 0; lev <= finest_level; ++lev ) {

         sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

         auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& rho_arr = ldata_p->state.const_array(mfi,DENSITY);
            auto const& sig_arr = sigma[lev]->array(mfi);
            amrex::ParallelFor(bx, [rho_arr,sig_arr,dummy_dt]
            AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
               sig_arr(i,j,k) = dummy_dt / rho_arr(i,j,k);
            });
         }
         scaleProj_RZ(lev,*sigma[lev]);
      }
   }

   // Set the velocity to the gravity field
   Vector<MultiFab> vel(finest_level+1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      vel[lev].define(grids[lev], dmap[lev], AMREX_SPACEDIM, nGhost, MFInfo(), *m_factory[lev]);
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
         vel[lev].setVal(m_gravity[idim],idim,1,1);
      }
      vel[lev].setBndry(0.0);
      setInflowBoundaryVel(vel[lev],lev,AmrNewTime);
      scaleProj_RZ(lev,vel[lev]);
   }

   // Done without divU in IAMR
   doNodalProject(GetVecOfPtrs(vel), GetVecOfPtrs(sigma), {}, {}, incremental, dummy_dt);
}

void PeleLM::velocityProjection(int is_initIter,
                                const TimeStamp &a_rhoTime,
                                const Real &a_dt)
{
   BL_PROFILE_VAR("PeleLM::velocityProjection()", velocityProjection);

   int nGhost = 0;
   int incremental = (is_initIter) ? 1 : 0;

   // Get sigma : scaled density inv. if not incompressible
   Vector<std::unique_ptr<MultiFab> > sigma(finest_level+1);
   if (! m_incompressible ) {
      Vector<std::unique_ptr<MultiFab> > rhoHalf(finest_level+1);
      rhoHalf = getDensityVect(a_rhoTime);
      for (int lev = 0; lev <= finest_level; ++lev ) {

         sigma[lev].reset(new MultiFab(grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

#ifdef AMREX_USE_OMP
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
#ifdef AMREX_USE_EB
         EB_set_covered(*sigma[lev],0.0);
#endif
         scaleProj_RZ(lev,*sigma[lev]);
      }
   }

   if (!incremental) {
      Vector<std::unique_ptr<MultiFab> > rhoHalf(finest_level+1);
      if (!m_incompressible) rhoHalf = getDensityVect(a_rhoTime);
      for (int lev = 0; lev <= finest_level; ++lev ) {

         auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
         auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
         for (MFIter mfi(ldataNew_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            Box const& bx = mfi.tilebox();
            auto const& vel_arr = ldataNew_p->state.array(mfi,VELX);
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
         MultiFab::Subtract(ldataNew_p->state,
                            ldataOld_p->state,
                            VELX,VELX,AMREX_SPACEDIM,0);
      }
   }

   // Get velocity
   Vector<std::unique_ptr<MultiFab> > vel;
   for (int lev = 0; lev <= finest_level; ++lev) {
      vel.push_back(std::make_unique<MultiFab> (m_leveldata_new[lev]->state,amrex::make_alias,VELX,AMREX_SPACEDIM));
#ifdef AMREX_USE_EB
      EB_set_covered(*vel[lev],0.0);
#endif
      vel[lev]->setBndry(0.0);
      if (!incremental) setInflowBoundaryVel(*vel[lev],lev,AmrNewTime);
      scaleProj_RZ(lev,*vel[lev]);
   }

   // To ensure integral of RHS is zero for closed chamber, get mean divU
   Real SbarOld = 0.0;
   Real SbarNew = 0.0;
   if (m_closed_chamber && !m_incompressible) {
      SbarNew = MFSum(GetVecOfConstPtrs(getDivUVect(AmrNewTime)),0);
      SbarNew /= m_uncoveredVol;        // Transform in Mean.
      if (incremental) {
         SbarOld = MFSum(GetVecOfConstPtrs(getDivUVect(AmrOldTime)),0);
         SbarOld /= m_uncoveredVol;        // Transform in Mean.
      }
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
            if (m_closed_chamber) {
               rhs_cc[lev].plus(-SbarNew,0,1);
            }
            rhs_cc[lev].mult(-1.0,0,1,ldata_p->divu.nGrow());
         } else {
            auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
            auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
            rhs_cc[lev].define(grids[lev],dmap[lev],1,ldataOld_p->divu.nGrow(), MFInfo(), *m_factory[lev]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
            for (MFIter mfi(rhs_cc[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
               const Box& gbx       = mfi.growntilebox();
               const auto& divu_o   = ldataOld_p->divu.const_array(mfi);
               const auto& divu_n   = ldataNew_p->divu.const_array(mfi);
               const auto& rhs      = rhs_cc[lev].array(mfi);
               amrex::ParallelFor(gbx, [divu_o,divu_n,rhs,SbarNew,SbarOld,is_closed_ch=m_closed_chamber]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
               {
                  rhs(i,j,k) = - (divu_n(i,j,k) - divu_o(i,j,k));
                  if (is_closed_ch) {
                     rhs(i,j,k) += SbarNew - SbarOld;               // subtract the mean, but rhs's already -
                  }
               });
            }
         }
#ifdef AMREX_USE_EB
         EB_set_covered(rhs_cc[lev],0.0);
#endif
         scaleProj_RZ(lev,rhs_cc[lev]);
      }
   }

   doNodalProject(GetVecOfPtrs(vel), GetVecOfPtrs(sigma), GetVecOfPtrs(rhs_cc), {}, incremental, a_dt);

   // If incremental
   // define back to be U^{np1} by adding U^{n}
   // and handles scaling if 2D-RZ
   if (incremental) {
      for (int lev = 0; lev <= finest_level; ++lev) {
         auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
         auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);
         unscaleProj_RZ(lev,*vel[lev]);         // Unscaling New vel before adding back old one
         MultiFab::Add(ldataNew_p->state,
                       ldataOld_p->state,
                       VELX,VELX,AMREX_SPACEDIM,0);
      }
   } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
         unscaleProj_RZ(lev,*vel[lev]);
      }
   }

   // In R-Z, AMReX-Hydro do an average down of r*vel.
   // Now that we have unscaled vel, need to do average down again
   // to have consistent vel across levels
   if (Geom(0).IsRZ()) {
      averageDownVelocity(AmrNewTime);
   }
}

void PeleLM::doNodalProject(const Vector<MultiFab*> &a_vel,
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
   std::unique_ptr<Hydro::NodalProjector> nodal_projector;

   if ( m_incompressible ) {
      Real constant_sigma = scaling_factor / m_rho;
      nodal_projector.reset(new Hydro::NodalProjector(a_vel, constant_sigma, Geom(0,finest_level), info));
   } else {
      if ( has_rhs ) {
         nodal_projector.reset(new Hydro::NodalProjector(a_vel, GetVecOfConstPtrs(a_sigma), Geom(0,finest_level),
                                                         info, rhs_cc, rhs_nd));
      } else {
         nodal_projector.reset(new Hydro::NodalProjector(a_vel, GetVecOfConstPtrs(a_sigma), Geom(0,finest_level), info));
      }
   }

   nodal_projector->setDomainBC(lobc, hibc);

#ifdef AMREX_USE_HYPRE
   nodal_projector->getMLMG().setHypreOptionsNamespace(m_hypre_namespace_nodal);
#endif

   // Solve
   nodal_projector->project(m_nodal_mg_rtol, m_nodal_mg_atol);

   auto phi = nodal_projector->getPhi();
   auto gphi = nodal_projector->getGradPhi();

   for(int lev = 0; lev <= finest_level; lev++) {

      auto ldata_p = getLevelDataPtr(lev,AmrNewTime);

#ifdef AMREX_USE_OMP
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

void
PeleLM::scaleProj_RZ(int a_lev,
                    MultiFab &a_mf)
{
#if AMREX_SPACEDIM == 2
    // Scale nodal projection cell-centered mfs by radius
    if (geom[a_lev].IsRZ()) {
        Box domain = geom[a_lev].Domain();
        auto BCRecVel = fetchBCRecArray(VELX,1);
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (BCRecVel[0].lo(idim) == BCType::ext_dir) {
                domain.growLo(idim,1);
            }
            if (BCRecVel[0].hi(idim) == BCType::ext_dir) {
                domain.growHi(idim,1);
            }
        }
        const Real dr     = geom[a_lev].CellSize()[0];
        auto const& mf_ma = a_mf.arrays();
        amrex::ParallelFor(a_mf, a_mf.nGrowVect(), [=,ncomp=a_mf.nComp()]
        AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            auto mf = mf_ma[box_no];
            if (domain.contains(i,j,k)) {
                for (int n = 0; n < ncomp; ++n) {
                   mf(i,j,k,n) *= (static_cast<Real>(i) + 0.5) * dr;
                }
            } else {
                for (int n = 0; n < ncomp; ++n) {
                   mf(i,j,k,n) = 0.0;
                }
            }
        });
    }
#endif
}

void
PeleLM::unscaleProj_RZ(int a_lev,
                       MultiFab &a_mf)
{
#if AMREX_SPACEDIM == 2
    // Unscale nodal projection cell-centered mfs by radius
    if (geom[a_lev].IsRZ()) {
        const Box& domain  = geom[a_lev].Domain();
        const Real dr      = geom[a_lev].CellSize()[0];
        auto const& mf_ma = a_mf.arrays();
        amrex::ParallelFor(a_mf, a_mf.nGrowVect(), [=,ncomp=a_mf.nComp()]
        AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
        {
            auto mf = mf_ma[box_no];
            if (domain.contains(i,j,k)) {
                for (int n = 0; n < ncomp; ++n) {
                   mf(i,j,k,n) /= (static_cast<Real>(i) + 0.5) * dr;
                }
            } else {
                for (int n = 0; n < ncomp; ++n) {
                   mf(i,j,k,n) = 0.0;
                }
            }
        });
    }
#endif
}
