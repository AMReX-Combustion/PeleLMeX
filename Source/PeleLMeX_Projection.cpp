#include <PeleLMeX.H>
#include <memory>

using namespace amrex;

void
PeleLM::initialProjection()
{
  BL_PROFILE("PeleLMeX::initialProjection()");

  if (m_verbose != 0) {
    Vector<Real> velMax(AMREX_SPACEDIM);
    velMax = MLNorm0(
      GetVecOfConstPtrs(getVelocityVect(AmrNewTime)), 0, AMREX_SPACEDIM);
    amrex::Print() << " Initial velocity projection: ";
    amrex::Print() << AMREX_D_TERM(
      "  U: " << velMax[0] <<, "  V: " << velMax[1] <<,
      "  W: " << velMax[2] <<) "\n";
  }

  constexpr Real dummy_dt = 1.0;
  constexpr int incremental = 0;
  constexpr int nGhost = 0;

  // Get sigma : density if not incompressible
  Vector<std::unique_ptr<MultiFab>> sigma;
  if (m_incompressible == 0) {
    sigma.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {

      sigma.emplace_back(
        std::make_unique<MultiFab>(
          grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

      auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
      auto state_ma = ldata_p->state.const_arrays();
      auto sigma_ma = sigma[lev]->arrays();
      amrex::ParallelFor(
        ldata_p->state, [state_ma, sigma_ma] AMREX_GPU_DEVICE(
                          int box_no, int i, int j, int k) noexcept {
          Array4<Real const> rho(state_ma[box_no], DENSITY);
          sigma_ma[box_no](i, j, k) = dummy_dt / rho(i, j, k);
        });
      Gpu::streamSynchronize();
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        scaleProj_RZ(lev, *sigma[lev]);
      }
#endif
    }
  }

  // Get velocity
  Vector<std::unique_ptr<MultiFab>> vel;
  vel.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    vel.emplace_back(
      std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, VELX, AMREX_SPACEDIM));
    vel[lev]->setBndry(0.0);
    setInflowBoundaryVel(*vel[lev], lev, AmrNewTime);
#if AMREX_SPACEDIM == 2
    if (geom[lev].IsRZ()) {
      scaleProj_RZ(lev, *vel[lev]);
    };
#endif
  }

  // Get RHS cc: - divU (- \int{divU})
  Real Sbar = 0.0;
  Vector<MultiFab> rhs_cc;
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
      MultiFab::Copy(
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
    Vector<Real> velMax = MLNorm0(
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

  constexpr Real dummy_dt = 1.0;
  constexpr int incremental = 0;
  constexpr int nGhost = 1;

  // Get sigma : density if not incompressible
  Vector<std::unique_ptr<MultiFab>> sigma;
  if (m_incompressible == 0) {
    sigma.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {

      sigma.emplace_back(
        std::make_unique<MultiFab>(
          grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

      auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
      auto state_ma = ldata_p->state.const_arrays();
      auto sigma_ma = sigma[lev]->arrays();
      amrex::ParallelFor(
        ldata_p->state, [state_ma, sigma_ma] AMREX_GPU_DEVICE(
                          int box_no, int i, int j, int k) noexcept {
          Array4<Real const> rho(state_ma[box_no], DENSITY);
          sigma_ma[box_no](i, j, k) = dummy_dt / rho(i, j, k);
        });
      Gpu::streamSynchronize();
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        scaleProj_RZ(lev, *sigma[lev]);
      }
#endif
    }
  }

  // Set the velocity to the gravity field
  Vector<MultiFab> vel;
  vel.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    vel.emplace_back(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGhost, MFInfo(), *m_factory[lev]);
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
  const int is_initIter, const TimeStamp a_rhoTime, const Real a_dt)
{
  BL_PROFILE("PeleLMeX::velocityProjection()");

  constexpr int nGhost = 0;
  const int incremental = (is_initIter) != 0 ? 1 : 0;

  // Get sigma : scaled density inv. if not incompressible
  Vector<std::unique_ptr<MultiFab>> sigma;
  if (m_incompressible == 0) {
    Vector<std::unique_ptr<MultiFab>> rhoHalf = getDensityVect(a_rhoTime);
    sigma.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {

      sigma.emplace_back(
        std::make_unique<MultiFab>(
          grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

      auto rhoHalf_ma = rhoHalf[lev]->const_arrays();
      auto sigma_ma = sigma[lev]->arrays();

      amrex::ParallelFor(
        *rhoHalf[lev], [rhoHalf_ma, sigma_ma, dt = a_dt] AMREX_GPU_DEVICE(
                         int box_no, int i, int j, int k) noexcept {
          sigma_ma[box_no](i, j, k) = dt / rhoHalf_ma[box_no](i, j, k);
        });

#ifdef AMREX_USE_EB
      EB_set_covered(*sigma[lev], 0.0);
#endif
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        scaleProj_RZ(lev, *sigma[lev]);
      }
#endif
    }
  }

  if (incremental == 0) {
    Vector<std::unique_ptr<MultiFab>> rhoHalf;
    if (m_incompressible == 0) {
      rhoHalf = getDensityVect(a_rhoTime);
    }
    rhoHalf.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {

      auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
      auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);

      auto state_old_ma = ldataOld_p->state.arrays();
      auto gp_new_ma = ldataNew_p->gp.const_arrays();
      if (m_incompressible == 0) {
        auto rho_ma = rhoHalf[lev]->const_arrays();
        amrex::ParallelFor(
          ldataNew_p->state,
          [state_old_ma, gp_new_ma, rho_ma, dt = a_dt] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            Array4<Real> vel(state_old_ma[box_no], VELX);
            const Real soverrho = dt / rho_ma[box_no](i, j, k);
            for (int n = 0; n < NUM_SPECIES; ++n) {
              vel(i, j, k, n) += gp_new_ma[box_no](i, j, k, n) * soverrho;
            }
          });
      } else {
        amrex::ParallelFor(
          ldataNew_p->state,
          [state_old_ma, gp_new_ma, rho = m_rho, dt = a_dt] AMREX_GPU_DEVICE(
            int box_no, int i, int j, int k) noexcept {
            Array4<Real> vel(state_old_ma[box_no], VELX);
            const Real soverrho = dt / rho;
            for (int n = 0; n < NUM_SPECIES; ++n) {
              vel(i, j, k, n) += gp_new_ma[box_no](i, j, k, n) * soverrho;
            }
          });
      }
    }
    Gpu::streamSynchronize();
  }

  // If incremental
  // define "vel" to be U^{np1*} - U^{n} rather than U^{np1*}
  if (incremental != 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
      auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
      MultiFab::Subtract(
        ldataNew_p->state, ldataOld_p->state, VELX, VELX, AMREX_SPACEDIM, 0);
    }
  }

  // Get velocity
  Vector<std::unique_ptr<MultiFab>> vel;
  vel.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    vel.emplace_back(
      std::make_unique<MultiFab>(
        m_leveldata_new[lev]->state, amrex::make_alias, VELX, AMREX_SPACEDIM));
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
  Real SbarOld = 0.0;
  Real SbarNew = 0.0;
  if ((m_closed_chamber != 0) && (m_incompressible == 0)) {
    SbarNew = MFSum(GetVecOfConstPtrs(getDivUVect(AmrNewTime)), 0);
    SbarNew /= m_uncoveredVol; // Transform in Mean.
    if (incremental != 0) {
      SbarOld = MFSum(GetVecOfConstPtrs(getDivUVect(AmrOldTime)), 0);
      SbarOld /= m_uncoveredVol; // Transform in Mean.
    }
  }

  // Get RHS cc
  Vector<MultiFab> rhs_cc;
  if ((m_incompressible == 0) && (m_has_divu != 0)) {
    rhs_cc.reserve(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      if (incremental == 0) {
        auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
        rhs_cc.emplace_back(
          grids[lev], dmap[lev], 1, ldata_p->divu.nGrow(), MFInfo(),
          *m_factory[lev]);
        MultiFab::Copy(
          rhs_cc[lev], ldata_p->divu, 0, 0, 1, ldata_p->divu.nGrow());
        if (m_closed_chamber != 0) {
          rhs_cc[lev].plus(-SbarNew, 0, 1);
        }
        rhs_cc[lev].mult(-1.0, 0, 1, ldata_p->divu.nGrow());
      } else {
        auto* ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
        auto* ldataNew_p = getLevelDataPtr(lev, AmrNewTime);
        rhs_cc[lev].define(
          grids[lev], dmap[lev], 1, ldataOld_p->divu.nGrow(), MFInfo(),
          *m_factory[lev]);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(rhs_cc[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
          const Box& gbx = mfi.growntilebox();
          const auto& divu_o = ldataOld_p->divu.const_array(mfi);
          const auto& divu_n = ldataNew_p->divu.const_array(mfi);
          const auto& rhs = rhs_cc[lev].array(mfi);
          amrex::ParallelFor(
            gbx,
            [divu_o, divu_n, rhs, SbarNew, SbarOld,
             is_closed_ch =
               m_closed_chamber] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              rhs(i, j, k) = -(divu_n(i, j, k) - divu_o(i, j, k));
              if (is_closed_ch != 0) {
                rhs(i, j, k) +=
                  SbarNew - SbarOld; // subtract the mean, but rhs's already -
              }
            });
        }
      }
#ifdef AMREX_USE_EB
      EB_set_covered(rhs_cc[lev], 0.0);
#endif
#if AMREX_SPACEDIM == 2
      if (geom[lev].IsRZ()) {
        scaleProj_RZ(lev, rhs_cc[lev]);
      }
#endif
    }
  }

  doNodalProject(
    GetVecOfPtrs(vel), GetVecOfPtrs(sigma), GetVecOfPtrs(rhs_cc), {},
    incremental, a_dt);

#if AMREX_SPACEDIM == 2
  for (int lev = 0; lev <= finest_level; ++lev) {
    // Unscaling New vel before adding back old one
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
      MultiFab::Add(
        ldataNew_p->state, ldataOld_p->state, VELX, VELX, AMREX_SPACEDIM, 0);
    }
  }

#ifdef AMREX_SPACEDIM == 2
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
  const Vector<MultiFab*>& a_vel,
  const Vector<MultiFab*>& a_sigma,
  const Vector<MultiFab*>& rhs_cc,
  const Vector<const MultiFab*>& rhs_nd,
  const int incremental,
  const Real scaling_factor)
{
  // Asserts
  AMREX_ASSERT(a_vel.size() == a_sigma.size());
  AMREX_ASSERT(rhs_cc.empty() || (a_vel.size() == rhs_cc.size()));
  AMREX_ASSERT(rhs_nd.empty() || (a_vel.size() == rhs_nd.size()));
  AMREX_ASSERT(a_vel[0]->nComp() == AMREX_SPACEDIM);

  LPInfo info;
  info.setMaxCoarseningLevel(m_nodal_mg_max_coarsening_level);

  // BCs
  std::array<LinOpBCType, AMREX_SPACEDIM> lobc;
  std::array<LinOpBCType, AMREX_SPACEDIM> hibc;
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (Geom(0).isPeriodic(idim)) {
      lobc[idim] = hibc[idim] = LinOpBCType::Periodic;
    } else {
      if (m_phys_bc.lo(idim) == amrex::PhysBCType::outflow) {
        lobc[idim] = LinOpBCType::Dirichlet;
      } else if (m_phys_bc.lo(idim) == amrex::PhysBCType::inflow) {
        lobc[idim] = LinOpBCType::inflow;
      } else {
        lobc[idim] = LinOpBCType::Neumann;
      }
      if (m_phys_bc.hi(idim) == amrex::PhysBCType::outflow) {
        hibc[idim] = LinOpBCType::Dirichlet;
      } else if (m_phys_bc.hi(idim) == amrex::PhysBCType::inflow) {
        hibc[idim] = LinOpBCType::inflow;
      } else {
        hibc[idim] = LinOpBCType::Neumann;
      }
    }
  }

  // Setup NodalProjector
  std::unique_ptr<Hydro::NodalProjector> nodal_projector;

  if (m_incompressible != 0) {
    const Real constant_sigma = scaling_factor / m_rho;
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
  nodal_projector->getMLMG().setHypreOptionsNamespace(m_hypre_namespace_nodal);
#endif

  // Solve
  nodal_projector->project(m_nodal_mg_rtol, m_nodal_mg_atol);

  auto phi = nodal_projector->getPhi();
  auto gphi = nodal_projector->getGradPhi();

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldata_p->gp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Box const& tbx = mfi.tilebox();
      Box const& nbx = mfi.nodaltilebox();
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
PeleLM::scaleProj_RZ( // NOLINT(readability-convert-member-functions-to-static)
		      const int a_lev,
		      MultiFab& a_mf)
{
  // Scale nodal projection cell-centered mfs by radius
  Box domain = geom[a_lev].Domain();
  auto BCRecVel = fetchBCRecArray(VELX, 1);
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (BCRecVel[0].lo(idim) == BCType::ext_dir) {
      domain.growLo(idim, 1);
    }
    if (BCRecVel[0].hi(idim) == BCType::ext_dir) {
      domain.growHi(idim, 1);
    }
  }
  const Real dr = geom[a_lev].CellSize()[0];
  auto const& mf_ma = a_mf.arrays();
  amrex::ParallelFor(
    a_mf, a_mf.nGrowVect(),
    [mf_ma, dr, domain, ncomp = a_mf.nComp()] AMREX_GPU_DEVICE(
      int box_no, int i, int j, int k) noexcept {
      auto mf = mf_ma[box_no];
      if (domain.contains(i, j, k)) {
        for (int n = 0; n < ncomp; ++n) {
          mf(i, j, k, n) *= (static_cast<Real>(i) + 0.5) * dr;
        }
      } else {
        for (int n = 0; n < ncomp; ++n) {
          mf(i, j, k, n) = 0.0;
        }
      }
    });
  Gpu::streamSynchronize();
}

void
PeleLM::
  unscaleProj_RZ( // NOLINT(readability-convert-member-functions-to-static)
    const int a_lev,
    MultiFab& a_mf)
{
  // Unscale nodal projection cell-centered mfs by radius
  const Box& domain = geom[a_lev].Domain();
  const Real dr = geom[a_lev].CellSize()[0];
  auto const& mf_ma = a_mf.arrays();
  amrex::ParallelFor(
    a_mf, a_mf.nGrowVect(),
    [mf_ma, dr, domain, ncomp = a_mf.nComp()] AMREX_GPU_DEVICE(
      int box_no, int i, int j, int k) noexcept {
      auto mf = mf_ma[box_no];
      if (domain.contains(i, j, k)) {
        for (int n = 0; n < ncomp; ++n) {
          mf(i, j, k, n) /= (static_cast<Real>(i) + 0.5) * dr;
        }
      } else {
        for (int n = 0; n < ncomp; ++n) {
          mf(i, j, k, n) = 0.0;
        }
      }
    });
  Gpu::streamSynchronize();
}
#endif
