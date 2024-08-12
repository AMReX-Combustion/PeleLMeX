#include <PeleLMeX.H>
#include <PeleLMeX_K.H>

using namespace amrex;

// Return velocity forces scaled by rhoInv
// including grapP term if add_gradP = 1
// including divTau if input Vector not empty
void
PeleLM::getVelForces(
  const TimeStamp& a_time,
  const Vector<MultiFab*>& a_divTau,
  const Vector<MultiFab*>& a_velForce,
  int nGrowForce,
  int add_gradP)
{
  BL_PROFILE("PeleLMeX::getVelForces()");
  int has_divTau = static_cast<int>(!a_divTau.empty());

  for (int lev = 0; lev <= finest_level; ++lev) {
    if (has_divTau != 0) {
      getVelForces(a_time, lev, a_divTau[lev], a_velForce[lev], add_gradP);
    } else {
      getVelForces(a_time, lev, nullptr, a_velForce[lev], add_gradP);
    }
  }

  // FillPatch forces
  if (nGrowForce > 0) {
    fillpatch_forces(m_cur_time, a_velForce, nGrowForce);
  }
}

void
PeleLM::getVelForces(
  const TimeStamp& a_time,
  int lev,
  MultiFab* a_divTau,
  MultiFab* a_velForce,
  int add_gradP)
{

  // Get level data
  // TODO: the 1 here bypass getting halftime vel and return oldtime vel
  auto* ldata_p = getLevelDataPtr(lev, a_time, 1);

  // Get gradp: if m_t_old < 0.0, we are during initialization -> only NewTime
  // data initialized at this point
  auto* ldataGP_p = (m_t_old[lev] < 0.0) ? getLevelDataPtr(lev, AmrNewTime)
                                         : getLevelDataPtr(lev, AmrOldTime);

  Real time = getTime(lev, a_time);

  int has_divTau = static_cast<int>(a_divTau != nullptr);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*a_velForce, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const auto& bx = mfi.tilebox();
    FArrayBox DummyFab(bx, 1);
    const auto& vel_arr = ldata_p->state.const_array(mfi, VELX);
    const auto& rho_arr = (m_incompressible) != 0
                            ? DummyFab.array()
                            : ldata_p->state.const_array(mfi, DENSITY);
    const auto& rhoY_arr = (m_incompressible) != 0
                             ? DummyFab.array()
                             : ldata_p->state.const_array(mfi, FIRSTSPEC);
    const auto& rhoh_arr = (m_incompressible) != 0
                             ? DummyFab.array()
                             : ldata_p->state.const_array(mfi, RHOH);
    const auto& temp_arr = (m_incompressible) != 0
                             ? DummyFab.array()
                             : ldata_p->state.const_array(mfi, TEMP);
    const auto& extmom_arr = m_extSource[lev]->const_array(mfi, VELX);
    const auto& extrho_arr = m_extSource[lev]->const_array(mfi, DENSITY);
    const auto& force_arr = a_velForce->array(mfi);

    // Get other forces (gravity, ...)
    getVelForces(
      lev, bx, time, force_arr, vel_arr, rho_arr, rhoY_arr, rhoh_arr, temp_arr,
      extmom_arr, extrho_arr);

#ifdef PELE_USE_EFIELD
    const auto& phiV_arr = ldata_p->state.const_array(mfi, PHIV);
    const auto& ne_arr = ldata_p->state.const_array(mfi, NE);
    addLorentzVelForces(lev, bx, time, force_arr, rhoY_arr, phiV_arr, ne_arr);
#endif

    // Add pressure gradient and viscous forces (if req.) and scale by density.
    int is_incomp = m_incompressible;
    Real incomp_rho_inv = 1.0 / m_rho;
    if ((add_gradP != 0) || (has_divTau != 0)) {
      const auto& gp_arr =
        (add_gradP) != 0 ? ldataGP_p->gp.const_array(mfi) : DummyFab.array();
      const auto& divTau_arr =
        (has_divTau) != 0 ? a_divTau->const_array(mfi) : DummyFab.array();
      amrex::ParallelFor(
        bx,
        [incomp_rho_inv, is_incomp, add_gradP, has_divTau, rho_arr, gp_arr,
         divTau_arr, force_arr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (is_incomp != 0) {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
              if (add_gradP != 0) {
                force_arr(i, j, k, idim) -= gp_arr(i, j, k, idim);
              }
              if (has_divTau != 0) {
                force_arr(i, j, k, idim) += divTau_arr(i, j, k, idim);
              }
              force_arr(i, j, k, idim) *= incomp_rho_inv;
            }
          } else {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
              if (add_gradP != 0) {
                force_arr(i, j, k, idim) -= gp_arr(i, j, k, idim);
              }
              if (has_divTau != 0) {
                force_arr(i, j, k, idim) += divTau_arr(i, j, k, idim);
              }
              force_arr(i, j, k, idim) /= rho_arr(i, j, k);
            }
          }
        });
    } else {
      amrex::ParallelFor(
        bx, [incomp_rho_inv, is_incomp, rho_arr,
             force_arr] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (is_incomp != 0) {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
              force_arr(i, j, k, idim) *= incomp_rho_inv;
            }
          } else {
            for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
              force_arr(i, j, k, idim) /= rho_arr(i, j, k);
            }
          }
        });
    }
  }
}

void
PeleLM::getVelForces(
  int lev,
  const Box& bx,
  const Real& a_time,
  Array4<Real> const& force,
  Array4<const Real> const& vel,
  Array4<const Real> const& rho,
  Array4<const Real> const& rhoY,
  Array4<const Real> const& rhoh,
  Array4<const Real> const& temp,
  Array4<const Real> const& extMom,
  Array4<const Real> const& extRho)
{
  const auto dx = geom[lev].CellSizeArray();

  // Get non-static info for the pseudo gravity forcing
  int pseudo_gravity = m_ctrl_pseudoGravity;
  const Real dV_control = m_ctrl_dV;

  int is_incomp = m_incompressible;
  Real rho_incomp = m_rho;

  amrex::ParallelFor(
    bx,
    [=, grav = m_gravity, gp0 = m_background_gp,
     ps_dir = m_ctrl_flameDir] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      makeVelForce(
        i, j, k, is_incomp, rho_incomp, pseudo_gravity, ps_dir, a_time, grav,
        gp0, dV_control, dx, vel, rho, rhoY, rhoh, temp, extMom, extRho, force);
    });
}

void
PeleLM::addSpark(const int lev, const TimeStamp& a_timestamp)
{

  const Real* probLo = geom[lev].ProbLo();
  auto const dx = geom[lev].CellSizeArray();
  for (int n = 0; n < m_n_sparks; n++) {
    IntVect spark_idx;
    Real time = getTime(lev, a_timestamp);
    if (
      time < m_spark_time[n] || time > m_spark_time[n] + m_spark_duration[n]) {
      if (m_spark_verbose > 1 && lev == 0 && a_timestamp == AmrOldTime) {
        Print() << m_spark[n] << " not active" << std::endl;
      }
      continue;
    }
    if (m_spark_verbose > 1 && lev == 0 && a_timestamp == AmrOldTime) {
      Print() << m_spark[n] << " active" << std::endl;
    }
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      spark_idx[d] = (int)((m_spark_location[n][d] - probLo[d]) / dx[d]);
    }
    Box domainBox = geom[lev].Domain();
    // just a check
    if (
      !domainBox.contains(spark_idx) && lev == 0 && a_timestamp == AmrOldTime) {
      Warning(m_spark[n] + " not in domain!");
      continue;
    }
    auto* ldata_p = getLevelDataPtr(lev, a_timestamp);
    auto eos = pele::physics::PhysicsType::eos();
    for (MFIter mfi(*(m_extSource[lev]), TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const Box& src_bx = mfi.growntilebox();
      auto const& rho_a = ldata_p->state.const_array(
        mfi, DENSITY); // need to convert temp to rhoH
      auto const& rhoY_a = ldata_p->state.const_array(mfi, FIRSTSPEC);
      auto const& temp_src_a = m_extSource[lev]->array(mfi, TEMP);
      auto const& rhoh_src_a = m_extSource[lev]->array(mfi, RHOH);
      amrex::ParallelFor(
        src_bx, [spark_duration = m_spark_duration[n], spark_idx,
                 spark_temp = m_spark_temp[n], spark_radius = m_spark_radius[n],
                 rho_a, rhoY_a, temp_src_a, rhoh_src_a, eos,
                 dx] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          Real dist_to_center = std::sqrt(AMREX_D_TERM(
            (i - spark_idx[0]) * (i - spark_idx[0]) * dx[0] * dx[0],
            +(j - spark_idx[1]) * (j - spark_idx[1]) * dx[1] * dx[1],
            +(k - spark_idx[2]) * (k - spark_idx[2]) * dx[2] * dx[2]));
          if (dist_to_center < spark_radius) {
            // probably doesn't actually
            // contribute to anything
            temp_src_a(i, j, k) = spark_temp / spark_duration;
            Real rhoh_src_loc = 0;
            Real rho = rho_a(i, j, k);
            Real Y[NUM_SPECIES];
            for (int ns = 0; ns < NUM_SPECIES; ns++) {
              Y[ns] = rhoY_a(i, j, k, ns) / rho;
            }
            eos.TY2H(spark_temp, Y, rhoh_src_loc);
            rhoh_src_loc *= rho * 1e-4 / spark_duration;
            rhoh_src_a(i, j, k) = rhoh_src_loc;
          }
        });
    }
  }
}
