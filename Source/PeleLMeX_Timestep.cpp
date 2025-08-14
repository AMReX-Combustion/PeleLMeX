#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

amrex::Real
PeleLM::computeDt(int is_init, const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::computeDt()");

  amrex::Real estdt = 1.0e200;

  //----------------------------------------------------------------
  // Store prev dt(s)
  m_prev_dt = m_dt;

  //----------------------------------------------------------------
  // Compute dt estimate from level data
  if (m_fixed_dt > 0.0) {
    estdt = m_fixed_dt;
  } else {
    if (((is_init != 0) || m_nstep == 0) && m_init_dt > 0.0) {
      estdt = m_init_dt;
    } else {
      amrex::Real dtconv = estConvectiveDt(a_time);
      estdt = amrex::min(estdt, dtconv);
      amrex::Real dtdivU = 1.0e200;
      if ((m_incompressible == 0) && (m_has_divu != 0)) {
        dtdivU = estDivUDt(a_time);
        estdt = amrex::min(estdt, dtdivU);
      }
#ifdef PELE_USE_PLASMA
      amrex::Real dtions = estEFIonsDt(a_time);
      estdt = amrex::min(estdt, dtions);
#endif
#ifdef PELE_USE_SPRAY
      amrex::Real dtspray = SprayEstDt();
      estdt = amrex::min(estdt, dtspray);
#endif
      if (m_verbose != 0) {
        amrex::Print() << " Est. time step - Conv: " << dtconv
                       << ", divu: " << dtdivU
#ifdef PELE_USE_PLASMA
                       << ", ions: " << dtions
#endif
#ifdef PELE_USE_SPRAY
                       << ", sprays: " << dtspray
#endif
                       << "\n";
      }
    }
  }

  //----------------------------------------------------------------
  // Limit dt
  if ((is_init != 0) || m_nstep == 0) {
    estdt *= m_dtshrink;
  } else {
    estdt = amrex::min(estdt, m_prev_dt * m_dtChangeMax);
    estdt = amrex::min(estdt, m_max_dt);
    // Shorten the dt to output plt file at exact req. time
    if (m_plot_per_exact > 0.0) {
      // Ensure ~O(dt) step by checking a little in advance
      amrex::Real timeToNextPlot =
        (std::floor(m_cur_time / m_plot_per_exact) + 1) * m_plot_per_exact -
        m_cur_time;
      if (2.0 * estdt > timeToNextPlot && timeToNextPlot > estdt) {
        estdt = amrex::Real(0.5) * timeToNextPlot;
      } else {
        if (timeToNextPlot > 1.e-12) {
          estdt = amrex::min(estdt, timeToNextPlot);
        }
      }
    }
    // If we're are getting close to the end of the simulation, shorten the dt
    // too
    if (m_stop_time >= 0.0) {
      // Ensure ~O(dt) last step by checking a little in advance
      amrex::Real timeLeft = (m_stop_time - m_cur_time);
      if (2.0 * estdt > timeLeft && timeLeft > estdt) {
        estdt = 0.5 * timeLeft;
      } else {
        estdt = amrex::min(estdt, timeLeft);
      }
    }
  }

  if (estdt < m_min_dt) {
    amrex::Print() << "\n";
    amrex::Print() << " ###################################### \n";
    amrex::Print() << " Estimated dt " << estdt << " is below allowed dt_min "
                   << m_min_dt << ": the simulation will stop ! \n";
    amrex::Print() << " ###################################### \n";
    amrex::Print() << "\n";
  }

  return estdt;
}

amrex::Real
PeleLM::estConvectiveDt(const TimeStamp& a_time)
{

  amrex::Real estdt = 1.0e200;
  constexpr amrex::Real small = 1.0e-8;

  for (int lev = 0; lev <= finest_level; ++lev) {

    amrex::Real estdt_lev = 1.0e200;

    //----------------------------------------------------------------
    // Get level data ptr
    auto* ldata_p = getLevelDataPtr(lev, a_time);

    auto const dx = geom[lev].CellSizeArray();

    //----------------------------------------------------------------
    // Get velocity forces
    int nGrow_force = 0;
    amrex::MultiFab velForces(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force, amrex::MFInfo(),
      Factory(lev));

    int add_gradP = 1;
    getVelForces(a_time, lev, nullptr, &velForces, add_gradP);

    //----------------------------------------------------------------
    // Get max forces
    amrex::Vector<amrex::Real> f_max(AMREX_SPACEDIM);
    f_max = velForces.norm0({AMREX_D_DECL(0, 1, 2)}, 0, true, true);

    // Get max velocity
    amrex::Vector<amrex::Real> u_max(AMREX_SPACEDIM);
    u_max =
      ldata_p->state.norm0({AMREX_D_DECL(VELX, VELY, VELZ)}, 0, true, true);

    //----------------------------------------------------------------
    // Est. min time step on lev
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (u_max[idim] > small) {
        estdt_lev = amrex::min(estdt_lev, dx[idim] / u_max[idim]);
      }
      if (f_max[idim] > small) {
        estdt_lev =
          amrex::min(estdt_lev, std::sqrt(2.0 * dx[idim] / f_max[idim]));
      }
    }

    //----------------------------------------------------------------
    // Set overall convective dt
    estdt = amrex::min(estdt, estdt_lev * m_cfl);
  }

  amrex::ParallelDescriptor::ReduceRealMin(estdt);

  return estdt;
}

amrex::Real
PeleLM::estDivUDt(const TimeStamp& a_time)
{

  amrex::Real estdt = 1.0e200;

  // Note: only methods 1 & 2 of PeleLM are available here
  AMREX_ASSERT(m_divu_checkFlag >= 0 && m_divu_checkFlag <= 2);

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, a_time);
    std::unique_ptr<amrex::MultiFab> density =
      std::make_unique<amrex::MultiFab>(
        ldata_p->state, amrex::make_alias, DENSITY, 1);

    auto dtfac = m_divu_dtFactor;
    auto rhoMin = m_divu_rhoMin;
    if (m_divu_checkFlag == 1) {
      amrex::Real divu_dt = amrex::ReduceMin(
        *density, ldata_p->divu, 0,
        [dtfac, rhoMin] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, amrex::Array4<amrex::Real const> const& rho,
          amrex::Array4<amrex::Real const> const& divu) -> amrex::Real {
          const auto lo = amrex::lbound(bx);
          const auto hi = amrex::ubound(bx);
          amrex::Real dt = 1.e37;
          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                amrex::Real dtcell =
                  est_divu_dt_1(i, j, k, dtfac, rhoMin, rho, divu);
                dt = amrex::min(dt, dtcell);
              }
            }
          }
          return dt;
        });
      estdt = amrex::min(divu_dt, estdt);
    } else if (m_divu_checkFlag == 2) {
      const auto& dxinv = geom[lev].InvCellSizeArray();
      std::unique_ptr<amrex::MultiFab> velo = std::make_unique<amrex::MultiFab>(
        ldata_p->state, amrex::make_alias, VELX, AMREX_SPACEDIM);
      amrex::Real divu_dt = amrex::ReduceMin(
        *density, ldata_p->divu, *velo, 0,
        [dtfac, rhoMin, dxinv] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, amrex::Array4<amrex::Real const> const& rho,
          amrex::Array4<amrex::Real const> const& vel,
          amrex::Array4<amrex::Real const> const& divu) -> amrex::Real {
          const auto lo = amrex::lbound(bx);
          const auto hi = amrex::ubound(bx);
          amrex::Real dt = 1.e37;
          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                amrex::Real dtcell =
                  est_divu_dt_2(i, j, k, dtfac, rhoMin, dxinv, rho, vel, divu);
                dt = amrex::min(dt, dtcell);
              }
            }
          }
          return dt;
        });
      estdt = amrex::min(divu_dt, estdt);
    }
  }

  amrex::ParallelDescriptor::ReduceRealMin(estdt);

  return estdt;
}

void
PeleLM::checkDt(const TimeStamp& a_time, const amrex::Real& a_dt)
{
  BL_PROFILE("PeleLMeX::checkDt()");

  if (m_fixed_dt > 0.0 || (m_divu_checkFlag == 0)) {
    return;
  }

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, a_time);

    const auto dxinv = geom[lev].InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(ldata_p->state, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& rho = ldata_p->state.const_array(mfi, DENSITY);
      auto const& vel = ldata_p->state.const_array(mfi, VELX);
      auto const& divu = ldata_p->divu.const_array(mfi);
      int divu_checkFlag = m_divu_checkFlag;
      auto dtfac = m_divu_dtFactor;
      auto rhoMin = m_divu_rhoMin;
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          check_divu_dt(
            i, j, k, divu_checkFlag, dtfac, rhoMin, dxinv, rho, vel, divu,
            a_dt);
        });
    }
  }
}
