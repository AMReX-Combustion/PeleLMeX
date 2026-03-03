#include <PeleLMeX.H>
#include <memory>

void
PeleLM::initialProjection()
{
  BL_PROFILE("PeleLMeX::initialProjection()");

  if (m_verbose != 0) {
    amrex::Vector<amrex::Real> velMax(AMREX_SPACEDIM);
    velMax = MLNorm0(
      GetVecOfConstPtrs(getVelocityVect(AmrNewTime)), 0, AMREX_SPACEDIM);
    amrex::Print() << " Initial velocity projection: ";
    amrex::Print() << AMREX_D_TERM(
      "  U: " << velMax[0] <<, "  V: " << velMax[1] <<,
      "  W: " << velMax[2] <<) "\n";
  }

  constexpr int incremental = 0;
  constexpr int nGhost = 0;
  constexpr amrex::Real dummy_dt = 1.0;

  // Get sigma : density if not incompressible
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> sigma(finest_level + 1);
  if (m_incompressible == 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      sigma[lev] = std::make_unique<amrex::MultiFab>(
        grids[lev], dmap[lev], 1, nGhost, amrex::MFInfo(), *m_factory[lev]);

      auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
      auto const& state_ma = ldata_p->state.const_arrays();
      auto const& sigma_ma = sigma[lev]->arrays();
      amrex::ParallelFor(
        ldata_p->state, [state_ma, sigma_ma] AMREX_GPU_DEVICE(
                          int box_no, int i, int j, int k) noexcept {
          amrex::Array4<amrex::Real const> rho(state_ma[box_no], DENSITY);
          sigma_ma[box_no](i, j, k) = dummy_dt / rho(i, j, k);
        });
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        amrex::Gpu::streamSynchronize();
        scaleProj_RZ(lev, *sigma[lev]);
      }
#endif
    }
    amrex::Gpu::streamSynchronize();
  }

  // Get velocity
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    vel[lev] = std::make_unique<amrex::MultiFab>(
      m_leveldata_new[lev]->state, amrex::make_alias, VELX, AMREX_SPACEDIM);
    vel[lev]->setBndry(0.0);
    setInflowBoundaryVel(*vel[lev], lev, AmrNewTime);
#if AMREX_SPACEDIM == 2
    if (geom[lev].IsRZ()) {
      scaleProj_RZ(lev, *vel[lev]);
    };
#endif
  }

  // Get RHS cc: - divU (- \int{divU})
  amrex::Real Sbar = 0.0;
  amrex::Vector<amrex::MultiFab> rhs_cc;
  if ((m_incompressible == 0) && (m_has_divu != 0)) {
    // Ensure integral of RHS is zero for closed chamber
    if (m_closed_chamber != 0) {
      Sbar = MFSum(GetVecOfConstPtrs(getDivUVect(AmrNewTime)), 0);
      Sbar /= m_uncoveredVol; // Transform in Mean.
    }
    rhs_cc.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      rhs_cc.emplace_back(
        grids[lev], dmap[lev], 1, m_leveldata_new[lev]->divu.nGrow());
      amrex::MultiFab::Copy(
        rhs_cc[lev], m_leveldata_new[lev]->divu, 0, 0, 1,
        m_leveldata_new[lev]->divu.nGrow());
      if (m_closed_chamber != 0) {
        rhs_cc[lev].plus(-Sbar, 0, 1);
      }
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        scaleProj_RZ(lev, rhs_cc[lev]);
      }
#endif
      rhs_cc[lev].mult(-1.0, 0, 1, rhs_cc[lev].nGrow());
    }
  }

  doNodalProject(
    GetVecOfPtrs(vel), GetVecOfPtrs(sigma), GetVecOfPtrs(rhs_cc), {},
    incremental, dummy_dt);

  // Set back press and gpress to zero and restore divu
  // and rescale velocity if 2D-RZ
  for (int lev = 0; lev <= finest_level; ++lev) {
    auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
    ldata_p->press.setVal(0.0);
    ldata_p->gp.setVal(0.0);
    if ((m_incompressible == 0) && (m_has_divu != 0)) {
      m_leveldata_new[lev]->divu.mult(-1.0, 0, 1, rhs_cc[lev].nGrow());
      // Restore divU integral
      if (m_closed_chamber != 0) {
        m_leveldata_new[lev]->divu.plus(Sbar, 0, 1);
      }
    }
#if AMREX_SPACEDIM == 2
    if (geom[lev].IsRZ()) {
      unscaleProj_RZ(lev, *vel[lev]);
    }
#endif
  }

  // In R-Z, AMReX-Hydro do an average down of r*vel.
  // Now that we have unscaled vel, need to do average down again
  // to have consistent vel across levels
  if (Geom(0).IsRZ()) {
    averageDownVelocity(AmrNewTime);
  }

  if (m_verbose != 0) {
    amrex::Vector<amrex::Real> velMax = MLNorm0(
      GetVecOfConstPtrs(getVelocityVect(AmrNewTime)), 0, AMREX_SPACEDIM);
    amrex::Print() << " >> After initial velocity projection: ";
    amrex::Print() << AMREX_D_TERM(
      "  U: " << velMax[0] <<, "  V: " << velMax[1] <<,
      "  W: " << velMax[2] <<) "\n";
  }
}

void
PeleLM::initialPressProjection()
{
  BL_PROFILE("PeleLMeX::initialPressProjection()");

  if (m_verbose != 0) {
    amrex::Print() << " Initial pressure projection \n";
  }

  constexpr amrex::Real dummy_dt = 1.0;
  constexpr int incremental = 0;
  constexpr int nGhost = 1;

  // Get sigma : density if not incompressible
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> sigma;
  if (m_incompressible == 0) {
    sigma.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {

      sigma.emplace_back(
        std::make_unique<amrex::MultiFab>(
          grids[lev], dmap[lev], 1, nGhost, amrex::MFInfo(), *m_factory[lev]));

      auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
      auto const& state_ma = ldata_p->state.const_arrays();
      auto const& sigma_ma = sigma[lev]->arrays();
      amrex::ParallelFor(
        ldata_p->state, [state_ma, sigma_ma] AMREX_GPU_DEVICE(
                          int box_no, int i, int j, int k) noexcept {
          amrex::Array4<amrex::Real const> rho(state_ma[box_no], DENSITY);
          sigma_ma[box_no](i, j, k) = dummy_dt / rho(i, j, k);
        });
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        amrex::Gpu::streamSynchronize();
        scaleProj_RZ(lev, *sigma[lev]);
      }
#endif
    }
    amrex::Gpu::streamSynchronize();
  }

  // Set the velocity to the gravity field
  amrex::Vector<amrex::MultiFab> vel;
  vel.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    vel.emplace_back(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGhost, amrex::MFInfo(),
      *m_factory[lev]);
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      vel[lev].setVal(m_gravity[idim], idim, 1, 1);
    }
    vel[lev].setBndry(0.0);
    setInflowBoundaryVel(vel[lev], lev, AmrNewTime);
#if AMREX_SPACEDIM == 2
    if (geom[lev].IsRZ()) {
      scaleProj_RZ(lev, vel[lev]);
    }
#endif
  }

  // Done without divU in IAMR
  doNodalProject(
    GetVecOfPtrs(vel), GetVecOfPtrs(sigma), {}, {}, incremental, dummy_dt);
}

void
PeleLM::velocityProjection(
  const int is_initIter, const TimeStamp a_rhoTime, const amrex::Real a_dt)
{
  BL_PROFILE("PeleLMeX::velocityProjection()");

  constexpr int nGhost = 0;
  const int incremental = (is_initIter) != 0 ? 1 : 0;

  // Get sigma : scaled density inv. if not incompressible
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> sigma(finest_level + 1);
  if (m_incompressible == 0) {
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhoHalf =
      getDensityVect(a_rhoTime);
    for (int lev = 0; lev <= finest_level; ++lev) {

      sigma[lev] = std::make_unique<amrex::MultiFab>(
        grids[lev], dmap[lev], 1, nGhost, amrex::MFInfo(), *m_factory[lev]);

      auto const& rhoHalf_ma = rhoHalf[lev]->const_arrays();
      auto const& sigma_ma = sigma[lev]->arrays();

      amrex::ParallelFor(
        *rhoHalf[lev], [rhoHalf_ma, sigma_ma, dt = a_dt] AMREX_GPU_DEVICE(
                         int box_no, int i, int j, int k) noexcept {
          sigma_ma[box_no](i, j, k) = dt / rhoHalf_ma[box_no](i, j, k);
        });
#ifdef AMREX_USE_EB
      amrex::Gpu::streamSynchronize();
      EB_set_covered(*sigma[lev], 0.0);
#endif
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        amrex::Gpu::streamSynchronize();
        scaleProj_RZ(lev, *sigma[lev]);
      }
#endif
    }
    amrex::Gpu::streamSynchronize();
  }

  if (incremental == 0) {
    if (m_incompressible == 0) {
      amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhoHalf =
        getDensityVect(a_rhoTime);
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
        auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
        auto const& state_new_ma = ldataNew_p->state.arrays();
        auto const& gp_old_ma = ldataOld_p->gp.const_arrays();
        auto const& rho_ma = rhoHalf[lev]->const_arrays();
        amrex::ParallelFor(
          ldataNew_p->state, amrex::IntVect(0), AMREX_SPACEDIM,
          [state_new_ma, gp_old_ma, rho_ma, dt = a_dt] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k, int n) noexcept {
            amrex::Array4<amrex::Real> vel(state_new_ma[box_no], VELX);
            const amrex::Real soverrho = dt / rho_ma[box_no](i, j, k);
            vel(i, j, k, n) += gp_old_ma[box_no](i, j, k, n) * soverrho;
          });
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
        auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
        auto const& state_new_ma = ldataNew_p->state.arrays();
        auto const& gp_old_ma = ldataOld_p->gp.const_arrays();
        const amrex::Real soverrho = m_dt / m_rho;
        amrex::ParallelFor(
          ldataNew_p->state, amrex::IntVect(0), AMREX_SPACEDIM,
          [state_new_ma, gp_old_ma, soverrho] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k, int n) noexcept {
            amrex::Array4<amrex::Real> vel(state_new_ma[box_no], VELX);
            vel(i, j, k, n) += gp_old_ma[box_no](i, j, k, n) * soverrho;
          });
      }
    }
    amrex::Gpu::streamSynchronize();
  }

  // If incremental
  // define "vel" to be U^{np1*} - U^{n} rather than U^{np1*}
  if (incremental != 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
      auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
      amrex::MultiFab::Subtract(
        ldataNew_p->state, ldataOld_p->state, VELX, VELX, AMREX_SPACEDIM, 0);
    }
  }

  // Get velocity
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> vel(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    vel[lev] = std::make_unique<amrex::MultiFab>(
      m_leveldata_new[lev]->state, amrex::make_alias, VELX, AMREX_SPACEDIM);
#ifdef AMREX_USE_EB
    EB_set_covered(*vel[lev], 0.0);
#endif
    vel[lev]->setBndry(0.0);
    if (incremental == 0) {
      setInflowBoundaryVel(*vel[lev], lev, AmrNewTime);
    }
#if AMREX_SPACEDIM == 2
    if (geom[lev].IsRZ()) {
      scaleProj_RZ(lev, *vel[lev]);
    }
#endif
  }

  // To ensure integral of RHS is zero for closed chamber, get mean divU
  amrex::Real SbarOld = 0.0;
  amrex::Real SbarNew = 0.0;
  if ((m_closed_chamber != 0) && (m_incompressible == 0)) {
    SbarNew = MFSum(GetVecOfConstPtrs(getDivUVect(AmrNewTime)), 0) /
              m_uncoveredVol; // Transform in Mean.
    if (incremental != 0) {
      SbarOld = MFSum(GetVecOfConstPtrs(getDivUVect(AmrOldTime)), 0) /
                m_uncoveredVol; // Transform in Mean.
    }
  }

  // Get RHS cc
  amrex::Vector<amrex::MultiFab> rhs_cc;
  if ((m_incompressible == 0) && (m_has_divu != 0)) {
    rhs_cc.reserve(finest_level + 1);
    if (incremental == 0) {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
        rhs_cc.emplace_back(
          grids[lev], dmap[lev], 1, ldata_p->divu.nGrow(), amrex::MFInfo(),
          *m_factory[lev]);
        amrex::MultiFab::Copy(
          rhs_cc[lev], ldata_p->divu, 0, 0, 1, ldata_p->divu.nGrow());
        if (m_closed_chamber != 0) {
          rhs_cc[lev].plus(-SbarNew, 0, 1);
        }
        rhs_cc[lev].mult(-1.0, 0, 1, ldata_p->divu.nGrow());
#ifdef AMREX_USE_EB
        EB_set_covered(rhs_cc[lev], 0.0);
#endif
#if AMREX_SPACEDIM == 2
        if (geom[lev].IsRZ()) {
          scaleProj_RZ(lev, rhs_cc[lev]);
        }
#endif
      }
    } else {
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
        auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
        rhs_cc.emplace_back(
          grids[lev], dmap[lev], 1, ldataOld_p->divu.nGrow(), amrex::MFInfo(),
          *m_factory[lev]);

        auto const& divu_o_ma = ldataOld_p->divu.const_arrays();
        auto const& divu_n_ma = ldataNew_p->divu.const_arrays();
        auto const& rhs_ma = rhs_cc[lev].arrays();

        amrex::ParallelFor(
          rhs_cc[lev], rhs_cc[lev].nGrowVect(),
          [divu_o_ma, divu_n_ma,
           rhs_ma] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
            rhs_ma[box_no](i, j, k) =
              -(divu_n_ma[box_no](i, j, k) - divu_o_ma[box_no](i, j, k));
          });
        if (m_closed_chamber != 0) {
          amrex::Gpu::streamSynchronize();
          rhs_cc[lev].plus(SbarNew - SbarOld, 0, 1, ldataOld_p->divu.nGrow());
        }
#ifdef AMREX_USE_EB
        amrex::Gpu::streamSynchronize();
        EB_set_covered(rhs_cc[lev], 0.0);
#endif
#if AMREX_SPACEDIM == 2
        if (geom[lev].IsRZ()) {
          amrex::Gpu::streamSynchronize();
          scaleProj_RZ(lev, rhs_cc[lev]);
        }
#endif
      }
      amrex::Gpu::streamSynchronize();
    }
  }

  doNodalProject(
    GetVecOfPtrs(vel), GetVecOfPtrs(sigma), GetVecOfPtrs(rhs_cc), {},
    incremental, a_dt);

#if AMREX_SPACEDIM == 2
  // Unscaling New vel before adding back old one
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (geom[lev].IsRZ()) {
      unscaleProj_RZ(lev, *vel[lev]);
    }
  }
#endif

  // If incremental
  // define back to be U^{np1} by adding U^{n}
  if (incremental != 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
      auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
      amrex::MultiFab::Add(
        ldataNew_p->state, ldataOld_p->state, VELX, VELX, AMREX_SPACEDIM, 0);
    }
  }

#if AMREX_SPACEDIM == 2
  // In R-Z, AMReX-Hydro do an average down of r*vel.
  // Now that we have unscaled vel, need to do average down again
  // to have consistent vel across levels
  if (Geom(0).IsRZ()) {
    averageDownVelocity(AmrNewTime);
  }
#endif
}

void
PeleLM::doNodalProject(
  const amrex::Vector<amrex::MultiFab*>& a_vel,
  const amrex::Vector<amrex::MultiFab*>& a_sigma,
  const amrex::Vector<amrex::MultiFab*>& rhs_cc,
  const amrex::Vector<const amrex::MultiFab*>& rhs_nd,
  const int incremental,
  const amrex::Real scaling_factor)
{
  // Asserts
  AMREX_ASSERT(a_vel.size() == a_sigma.size());
  AMREX_ASSERT(rhs_cc.empty() || (a_vel.size() == rhs_cc.size()));
  AMREX_ASSERT(rhs_nd.empty() || (a_vel.size() == rhs_nd.size()));
  AMREX_ASSERT(a_vel[0]->nComp() == AMREX_SPACEDIM);

  amrex::LPInfo info;
  info.setMaxCoarseningLevel(m_nodal_mg_max_coarsening_level);

  // BCs
  std::array<amrex::LinOpBCType, AMREX_SPACEDIM> lobc;
  std::array<amrex::LinOpBCType, AMREX_SPACEDIM> hibc;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (Geom(0).isPeriodic(idim)) {
      lobc[idim] = hibc[idim] = amrex::LinOpBCType::Periodic;
    } else {
      if (m_phys_bc.lo(idim) == amrex::PhysBCType::outflow) {
        lobc[idim] = amrex::LinOpBCType::Dirichlet;
      } else if (m_phys_bc.lo(idim) == amrex::PhysBCType::inflow) {
        lobc[idim] = amrex::LinOpBCType::inflow;
      } else {
        lobc[idim] = amrex::LinOpBCType::Neumann;
      }
      if (m_phys_bc.hi(idim) == amrex::PhysBCType::outflow) {
        hibc[idim] = amrex::LinOpBCType::Dirichlet;
      } else if (m_phys_bc.hi(idim) == amrex::PhysBCType::inflow) {
        hibc[idim] = amrex::LinOpBCType::inflow;
      } else {
        hibc[idim] = amrex::LinOpBCType::Neumann;
      }
    }
  }

  // Setup NodalProjector
  std::unique_ptr<Hydro::NodalProjector> nodal_projector;

  if (m_incompressible != 0) {
    const amrex::Real constant_sigma = scaling_factor / m_rho;
    nodal_projector = std::make_unique<Hydro::NodalProjector>(
      a_vel, constant_sigma, Geom(0, finest_level), info);
  } else {
    if (!rhs_cc.empty()) {
      nodal_projector = std::make_unique<Hydro::NodalProjector>(
        a_vel, GetVecOfConstPtrs(a_sigma), Geom(0, finest_level), info, rhs_cc,
        rhs_nd);
    } else {
      nodal_projector = std::make_unique<Hydro::NodalProjector>(
        a_vel, GetVecOfConstPtrs(a_sigma), Geom(0, finest_level), info);
    }
  }

  nodal_projector->setDomainBC(lobc, hibc);

#ifdef AMREX_USE_EB
  if (m_useEBinflow != 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      nodal_projector->getLinOp().setEBInflowVelocity(
        lev, *getEBState(lev, VELX, AMREX_SPACEDIM, AmrNewTime));
    }
  }
#endif

#ifdef AMREX_USE_HYPRE
  if (
    nodal_projector->getMLMG().getBottomSolver() ==
    amrex::MLMG::BottomSolver::hypre) {
    nodal_projector->getMLMG().setHypreOptionsNamespace(
      m_hypre_namespace_nodal);
  }
#endif

  // Solve
  if (!m_mlmg_fail_plt_residuals) {
    nodal_projector->project(m_nodal_mg_rtol, m_nodal_mg_atol);
  } else {
    nodal_projector->getMLMG().setThrowException(true);
    nodal_projector->getMLMG().setConvergenceNormType(
      amrex::MLMGNormType::bnorm);
    try {
      nodal_projector->project(m_nodal_mg_rtol, m_nodal_mg_atol);
    } catch (const std::exception& e) {
      amrex::Print() << "\n";
      amrex::Print() << "  *** Nodal projection MLMG solve failed! ***\n";
      amrex::Print() << "  Error: " << e.what() << "\n";
      amrex::Print()
        << "  Dumping nodal projection residuals for debugging...\n";

      auto phi = nodal_projector->getPhi();
      auto rhs = nodal_projector->getRHSConst();

      WriteMLMGResidual(
        nodal_projector->getMLMG(), phi, rhs, "nodal_projection", m_nstep);

      amrex::Abort("MLMG solve for nodal_projection failed");
    }
  }

  auto phi = nodal_projector->getPhi();
  auto gphi = nodal_projector->getGradPhi();

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(ldata_p->gp, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      amrex::Box const& tbx = mfi.tilebox();
      amrex::Box const& nbx = mfi.nodaltilebox();
      auto const& p_lev_arr = ldata_p->press.array(mfi);
      auto const& gp_lev_arr = ldata_p->gp.array(mfi);
      auto const& p_proj_arr = phi[lev]->const_array(mfi);
      auto const& gp_proj_arr = gphi[lev]->const_array(mfi);
      if (incremental != 0) {
        amrex::ParallelFor(
          tbx, AMREX_SPACEDIM,
          [gp_lev_arr,
           gp_proj_arr] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            gp_lev_arr(i, j, k, n) += gp_proj_arr(i, j, k, n);
          });
        amrex::ParallelFor(
          nbx, [p_lev_arr,
                p_proj_arr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            p_lev_arr(i, j, k) += p_proj_arr(i, j, k);
          });
      } else {
        amrex::ParallelFor(
          tbx, AMREX_SPACEDIM,
          [gp_lev_arr,
           gp_proj_arr] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            gp_lev_arr(i, j, k, n) = gp_proj_arr(i, j, k, n);
          });
        amrex::ParallelFor(
          nbx, [p_lev_arr,
                p_proj_arr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            p_lev_arr(i, j, k) = p_proj_arr(i, j, k);
          });
      }
    }
  }

  // Average down grad P
  for (int lev = finest_level - 1; lev >= 0; --lev) {
    auto* ldataFine_p = getLevelDataPtr(lev + 1, AmrNewTime);
    auto* ldataCrse_p = getLevelDataPtr(lev, AmrNewTime);
#ifdef AMREX_USE_EB
    amrex::EB_average_down(
      ldataFine_p->gp, ldataCrse_p->gp, 0, AMREX_SPACEDIM, refRatio(lev));
#else
    amrex::average_down(
      ldataFine_p->gp, ldataCrse_p->gp, 0, AMREX_SPACEDIM, refRatio(lev));
#endif
  }
}

#if AMREX_SPACEDIM == 2
void
PeleLM::scaleProj_RZ(const int a_lev, amrex::MultiFab& a_mf)
{
  // Scale nodal projection cell-centered mfs by radius
  amrex::Box domain = geom[a_lev].Domain();
  auto BCRecVel = fetchBCRecArray(VELX, 1);
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (BCRecVel[0].lo(idim) == amrex::BCType::ext_dir) {
      domain.growLo(idim, 1);
    }
    if (BCRecVel[0].hi(idim) == amrex::BCType::ext_dir) {
      domain.growHi(idim, 1);
    }
  }
  const amrex::Real dr = geom[a_lev].CellSize()[0];
  auto const& mf_ma = a_mf.arrays();
  amrex::ParallelFor(
    a_mf, a_mf.nGrowVect(), a_mf.nComp(),
    [mf_ma, dr,
     domain] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
      auto mf = mf_ma[box_no];
      if (domain.contains(i, j, k)) {
        mf(i, j, k, n) *= (static_cast<amrex::Real>(i) + 0.5) * dr;
      } else {
        mf(i, j, k, n) = 0.0;
      }
    });
  amrex::Gpu::streamSynchronize();
}

void
PeleLM::unscaleProj_RZ(const int a_lev, amrex::MultiFab& a_mf)
{
  // Unscale nodal projection cell-centered mfs by radius
  const amrex::Box& domain = geom[a_lev].Domain();
  const amrex::Real dr = geom[a_lev].CellSize()[0];
  auto const& mf_ma = a_mf.arrays();
  amrex::ParallelFor(
    a_mf, a_mf.nGrowVect(), a_mf.nComp(),
    [mf_ma, dr,
     domain] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
      auto mf = mf_ma[box_no];
      if (domain.contains(i, j, k)) {
        mf(i, j, k, n) /= (static_cast<amrex::Real>(i) + 0.5) * dr;
      } else {
        mf(i, j, k, n) = 0.0;
      }
    });
  amrex::Gpu::streamSynchronize();
}
#endif
