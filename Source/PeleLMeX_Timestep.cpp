#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

using namespace amrex;

Real
PeleLM::computeDt(int is_init, const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::computeDt()");

  Real estdt = 1.0e200;

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
      Real dtconv = estConvectiveDt(a_time);
      estdt = std::min(estdt, dtconv);
      Real dtdivU = 1.0e200;
      if ((m_incompressible == 0) && (m_has_divu != 0)) {
        dtdivU = estDivUDt(a_time);
        estdt = std::min(estdt, dtdivU);
      }
#ifdef PELE_USE_PLASMA
      Real dtions = estEFIonsDt(a_time);
      estdt = std::min(estdt, dtions);
#endif
#ifdef PELE_USE_SPRAY
      Real dtspray = SprayEstDt();
      estdt = std::min(estdt, dtspray);
#endif
      if (m_verbose != 0) {
        Print() << " Est. time step - Conv: " << dtconv << ", divu: " << dtdivU
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
    estdt = std::min(estdt, m_prev_dt * m_dtChangeMax);
    estdt = std::min(estdt, m_max_dt);
    // Shorten the dt to output plt file at exact req. time
    if (m_plot_per_exact > 0.0) {
      // Ensure ~O(dt) step by checking a little in advance
      Real timeToNextPlot =
        (std::floor(m_cur_time / m_plot_per_exact) + 1) * m_plot_per_exact -
        m_cur_time;
      if (2.0 * estdt > timeToNextPlot && timeToNextPlot > estdt) {
        estdt = Real(0.5) * timeToNextPlot;
      } else {
        if (timeToNextPlot > 1.e-12) {
          estdt = std::min(estdt, timeToNextPlot);
        }
      }
    }
    // If we're are getting close to the end of the simulation, shorten the dt
    // too
    if (m_stop_time >= 0.0) {
      // Ensure ~O(dt) last step by checking a little in advance
      Real timeLeft = (m_stop_time - m_cur_time);
      if (2.0 * estdt > timeLeft && timeLeft > estdt) {
        estdt = 0.5 * timeLeft;
      } else {
        estdt = std::min(estdt, timeLeft);
      }
    }
  }

  if (estdt < m_min_dt) {
    Print() << "\n";
    Print() << " ###################################### \n";
    Print() << " Estimated dt " << estdt << " is below allowed dt_min "
            << m_min_dt << ": the simulation will stop ! \n";
    Print() << " ###################################### \n";
    Print() << "\n";
  }

  return estdt;
}

Real
PeleLM::estConvectiveDt(const TimeStamp& a_time)
{

  Real estdt = 1.0e200;
  constexpr Real small = 1.0e-8;

  for (int lev = 0; lev <= finest_level; ++lev) {

    Real estdt_lev = 1.0e200;

    //----------------------------------------------------------------
    // Get level data ptr
    auto* ldata_p = getLevelDataPtr(lev, a_time);

    auto const dx = geom[lev].CellSizeArray();

    //----------------------------------------------------------------
    // Get velocity forces
    int nGrow_force = 0;
    MultiFab velForces(
      grids[lev], dmap[lev], AMREX_SPACEDIM, nGrow_force, MFInfo(),
      Factory(lev));

    int add_gradP = 1;
    getVelForces(a_time, lev, nullptr, &velForces, add_gradP);

    //----------------------------------------------------------------
    // Get max forces
    Vector<Real> f_max(AMREX_SPACEDIM);
    f_max = velForces.norm0({AMREX_D_DECL(0, 1, 2)}, 0, true, true);

    // Get max velocity
    Vector<Real> u_max(AMREX_SPACEDIM);
    u_max =
      ldata_p->state.norm0({AMREX_D_DECL(VELX, VELY, VELZ)}, 0, true, true);

    //----------------------------------------------------------------
    // Est. min time step on lev
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (u_max[idim] > small) {
        estdt_lev = std::min(estdt_lev, dx[idim] / u_max[idim]);
      }
      if (f_max[idim] > small) {
        estdt_lev =
          std::min(estdt_lev, std::sqrt(2.0 * dx[idim] / f_max[idim]));
      }
    }

    //----------------------------------------------------------------
    // Set overall convective dt
    estdt = std::min(estdt, estdt_lev * m_cfl);
  }

  ParallelDescriptor::ReduceRealMin(estdt);

  return estdt;
}

Real
PeleLM::estDivUDt(const TimeStamp& a_time)
{

  Real estdt = 1.0e200;

  // Note: only methods 1 & 2 of PeleLM are available here
  AMREX_ASSERT(m_divu_checkFlag >= 0 && m_divu_checkFlag <= 2);

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, a_time);
    std::unique_ptr<MultiFab> density =
      std::make_unique<MultiFab>(ldata_p->state, amrex::make_alias, DENSITY, 1);

    auto dtfac = m_divu_dtFactor;
    auto rhoMin = m_divu_rhoMin;
    if (m_divu_checkFlag == 1) {
      Real divu_dt = amrex::ReduceMin(
        *density, ldata_p->divu, 0,
        [dtfac, rhoMin] AMREX_GPU_HOST_DEVICE(
          Box const& bx, Array4<Real const> const& rho,
          Array4<Real const> const& divu) -> Real {
          using namespace amrex::literals;
          const auto lo = amrex::lbound(bx);
          const auto hi = amrex::ubound(bx);
          amrex::Real dt = 1.e37_rt;
          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                Real dtcell = est_divu_dt_1(i, j, k, dtfac, rhoMin, rho, divu);
                dt = amrex::min(dt, dtcell);
              }
            }
          }
          return dt;
        });
      estdt = std::min(divu_dt, estdt);
    } else if (m_divu_checkFlag == 2) {
      const auto& dxinv = geom[lev].InvCellSizeArray();
      std::unique_ptr<MultiFab> velo = std::make_unique<MultiFab>(
        ldata_p->state, amrex::make_alias, VELX, AMREX_SPACEDIM);
      Real divu_dt = amrex::ReduceMin(
        *density, ldata_p->divu, *velo, 0,
        [dtfac, rhoMin, dxinv] AMREX_GPU_HOST_DEVICE(
          Box const& bx, Array4<Real const> const& rho,
          Array4<Real const> const& vel,
          Array4<Real const> const& divu) -> Real {
          using namespace amrex::literals;
          const auto lo = amrex::lbound(bx);
          const auto hi = amrex::ubound(bx);
          amrex::Real dt = 1.e37_rt;
          for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
              for (int i = lo.x; i <= hi.x; ++i) {
                Real dtcell =
                  est_divu_dt_2(i, j, k, dtfac, rhoMin, dxinv, rho, vel, divu);
                dt = amrex::min(dt, dtcell);
              }
            }
          }
          return dt;
        });
      estdt = std::min(divu_dt, estdt);
    }
  }

  ParallelDescriptor::ReduceRealMin(estdt);

  return estdt;
}

void
PeleLM::checkDt(const TimeStamp& a_time, const Real& a_dt)
{
  BL_PROFILE("PeleLMeX::checkDt()");

  if (m_fixed_dt > 0.0 || (m_divu_checkFlag == 0)) {
    return;
  }

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, a_time);

    const auto dxinv = geom[lev].InvCellSizeArray();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldata_p->state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& rho = ldata_p->state.const_array(mfi, DENSITY);
      auto const& vel = ldata_p->state.const_array(mfi, VELX);
      auto const& divu = ldata_p->divu.const_array(mfi);
      int divu_checkFlag = m_divu_checkFlag;
      auto dtfac = m_divu_dtFactor;
      auto rhoMin = m_divu_rhoMin;
      amrex::ParallelFor(
        bx, [rho, vel, divu, divu_checkFlag, dtfac, rhoMin, dxinv,
             a_dt] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          check_divu_dt(
            i, j, k, divu_checkFlag, dtfac, rhoMin, dxinv, rho, vel, divu,
            a_dt);
        });
    }
  }
}
