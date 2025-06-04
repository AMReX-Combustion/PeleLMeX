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

#ifdef PELE_USE_PLASMA
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
PeleLM::addSpark(const TimeStamp& a_timestamp)
{
  for (int lev = 0; lev <= finest_level; lev++) {
    for (int n = 0; n < m_n_sparks; n++) {
      // Do the checks first
      Real time = getTime(lev, a_timestamp);
      bool verb = m_spark_verbose > 1 && lev == 0;
      if (
        time < m_spark_time[n] ||
        time > m_spark_time[n] + m_spark_duration[n]) {
        if (verb) {
          Print() << m_spark[n] << " not active" << std::endl;
        }
        continue;
      }
      const Real* probLo = geom[lev].ProbLo();
      auto const dx = geom[lev].CellSizeArray();
      IntVect spark_idx;
      for (int d = 0; d < AMREX_SPACEDIM; d++) {
        spark_idx[d] = (int)((m_spark_location[n][d] - probLo[d]) / dx[d]);
      }
      Box domainBox = geom[lev].Domain();
      // just a check
      if (!domainBox.contains(spark_idx)) {
        Warning(m_spark[n] + " not in domain!");
        continue;
      }
      if (verb) {
        Print() << m_spark[n] << " active" << std::endl;
      }

      auto statema = getLevelDataPtr(lev, a_timestamp)->state.const_arrays();
      auto extma = m_extSource[lev]->arrays();
      auto const* leosparm = eos_parms.device_parm();

      amrex::ParallelFor(
        *m_extSource[lev],
        [=, spark_duration = m_spark_duration[n], spark_temp = m_spark_temp[n],
         eosparm = leosparm,
         spark_radius = m_spark_radius
           [n]] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
          auto eos = pele::physics::PhysicsType::eos(eosparm);
          Real dist_to_center = std::sqrt(AMREX_D_TERM(
            (i - spark_idx[0]) * (i - spark_idx[0]) * dx[0] * dx[0],
            +(j - spark_idx[1]) * (j - spark_idx[1]) * dx[1] * dx[1],
            +(k - spark_idx[2]) * (k - spark_idx[2]) * dx[2] * dx[2]));
          if (dist_to_center < spark_radius) {
            Real rhoh_src_loc = 0;
            Real rho = statema[box_no](i, j, k, DENSITY);
            Real Y[NUM_SPECIES];
            for (int ns = 0; ns < NUM_SPECIES; ns++) {
              Y[ns] = statema[box_no](i, j, k, FIRSTSPEC + ns) / rho;
            }
            eos.TY2H(spark_temp, Y, rhoh_src_loc);
            rhoh_src_loc *= rho * 1e-4 / spark_duration;
            extma[box_no](i, j, k, RHOH) = rhoh_src_loc;
          }
        });
      Gpu::streamSynchronize();
    }
  }
}

// Manifold model - dissipation rate sources for variances
void
PeleLM::addScalarVarianceSources(const TimeStamp& a_timestamp)
{
  BL_PROFILE("PeleLM::addScalarVarianceSources");
  // no scalar dissipation sources if not using a manifold model
#ifndef USE_MANIFOLD_EOS
  amrex::ignore_unused(a_timestamp);
#else

  auto const& leosparm = eos_parms.host_parm();

  // Determine if we have any scalar variances that need sources added
  // Could be moved elsewhere to not do every timestep
  int nvariances = 0;
  int var_of_scalar = -1;
  for (int n = 0; n < MANIFOLD_DIM; ++n) {
    if (leosparm.is_variance_of[n] >= 0) {
      if (!m_do_les) {
        amrex::Abort("PeleLM::addScalarVarianceSources(): cannot add a "
                     "scalar dissipation without an active LES model");
      }
      nvariances += 1;
      var_of_scalar = FIRSTSPEC + leosparm.is_variance_of[n];
    }
  }

  if (nvariances > 1) {
    amrex::Abort(
      "PeleLM::addScalarVarianceSources(): currently we only support "
      "manifold models with 0 or 1 variances");
  } else if (nvariances > 0) {

    // Compute scalar gradients (no need to average down here)
    int do_avgDown = 0;
    auto bcRecScalar = fetchBCRecArray(var_of_scalar, 1);
    int nGrow = 0; // No need for ghost face on fluxes
    Vector<Array<MultiFab, AMREX_SPACEDIM>> grad_fc(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      const auto& ba = grids[lev];
      const auto& factory = Factory(lev);
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        grad_fc[lev][idim].define(
          amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
          nGrow, MFInfo(), factory);
        grad_fc[lev][idim].setVal(0.0); // Required?
      }
    }
    getDiffusionOp()->computeGradient(
      GetVecOfArrOfPtrs(grad_fc), {},
      GetVecOfConstPtrs(getStateVect(a_timestamp)), {}, bcRecScalar[0],
      do_avgDown, var_of_scalar);

    // Add in Production and Dissipation source terms for subfilter variances
    for (int lev = 0; lev <= finest_level; lev++) {

      auto* ldata_p = getLevelDataPtr(lev, a_timestamp);

      // Require the turbulent viscosity to be pre-computed
      // it always is stored at AmrOldTime, so we just use that
      // We need cell-centered mu_t but have it at faces
      // The simple interpolation below probably isn't valid for EB
#ifdef AMREX_USE_EB
      amrex::Abort(
        "PeleLM::addScalarVarianceSources(): this is not supported with EB");
#endif

      for (int n = 0; n < MANIFOLD_DIM; ++n) {
        if (leosparm.is_variance_of[n] >= 0) {

          constexpr amrex::Real fact = 0.5 / AMREX_SPACEDIM;
          const amrex::Real C_chi = m_les_c_chi;
          const amrex::Real ScInv = m_Schmidt_inv;

          AMREX_D_TERM(auto const& mut_arr_x =
                         m_leveldata_old[lev]->visc_turb_fc[0].const_arrays();
                       , auto const& mut_arr_y =
                           m_leveldata_old[lev]->visc_turb_fc[1].const_arrays();
                       ,
                       auto const& mut_arr_z =
                         m_leveldata_old[lev]->visc_turb_fc[2].const_arrays();)
          AMREX_D_TERM(auto const& gx = grad_fc[lev][0].const_arrays();
                       , auto const& gy = grad_fc[lev][1].const_arrays();
                       , auto const& gz = grad_fc[lev][2].const_arrays();)
          auto extma = m_extSource[lev]->arrays();
          auto statema = ldata_p->state.const_arrays();

          // l_scale will also need modification for EB
          const amrex::Real vol = AMREX_D_TERM(
            geom[lev].CellSize(0), *geom[lev].CellSize(1),
            *geom[lev].CellSize(2));
          const amrex::Real l_scale =
            (AMREX_SPACEDIM == 2) ? std::sqrt(vol) : std::cbrt(vol);
          const amrex::Real inv_l_scale2 = 1.0 / (l_scale * l_scale);

          amrex::ParallelFor(
            *m_extSource[lev],
            [=] AMREX_GPU_DEVICE(int bx, int i, int j, int k) noexcept {
              // Subfilter Scalar Dissipation: Linear Relaxation model
              // rho chi_sgs = C_chi * mu_t / Delta^2 * Variance
              amrex::Real mu_t =
                fact * (AMREX_D_TERM(
                         mut_arr_x[bx](i, j, k) + mut_arr_x[bx](i + 1, j, k),
                         +mut_arr_y[bx](i, j, k) + mut_arr_y[bx](i, j + 1, k),
                         +mut_arr_z[bx](i, j, k) + mut_arr_z[bx](i, j, k + 1)));

              extma[bx](i, j, k, FIRSTSPEC + n) -=
                C_chi * mu_t * inv_l_scale2 *
                statema[bx](i, j, k, FIRSTSPEC + n);

              // Production term (w/ Smagorinsky closure for turbulent flux)
              // -2 (rho <u_j C> - rho <u_j><C>) d<C>/dx_j
              // = 2 *mu_t/Sc_t * d<C>/dx_j * d<C>/dx_j
              amrex::Real mu_grad2 =
                fact *
                (AMREX_D_TERM(
                  mut_arr_x[bx](i, j, k) * gx[bx](i, j, k) * gx[bx](i, j, k) +
                    mut_arr_x[bx](i + 1, j, k) * gx[bx](i + 1, j, k) *
                      gx[bx](i + 1, j, k),
                  +mut_arr_y[bx](i, j, k) * gy[bx](i, j, k) * gy[bx](i, j, k) +
                    mut_arr_y[bx](i, j + 1, k) * gy[bx](i, j + 1, k) *
                      gy[bx](i, j + 1, k),
                  +mut_arr_z[bx](i, j, k) * gz[bx](i, j, k) * gz[bx](i, j, k) +
                    mut_arr_z[bx](i, j, k + 1) * gz[bx](i, j, k + 1) *
                      gz[bx](i, j, k + 1)));

              extma[bx](i, j, k, FIRSTSPEC + n) += 2.0 * ScInv * mu_grad2;
            });
        }
      }
      Gpu::streamSynchronize();
    }
  }

#endif
}

// Calculate additional external sources (soot, radiation, user defined, etc.)
void
PeleLM::getExternalSources(
  int is_initIter,
  const PeleLM::TimeStamp& a_timestamp_old,
  const PeleLM::TimeStamp& a_timestamp_new)
{
  amrex::ignore_unused(is_initIter);

  if (m_n_sparks > 0) {
    addSpark(a_timestamp_old);
  }

#ifdef PELE_USE_SPRAY
  if (is_initIter == 0) {
    SprayMKD(m_cur_time, m_dt);
  }
#endif
#ifdef PELE_USE_SOOT
  if (do_soot_solve) {
    computeSootSource(a_timestamp_old, m_dt);
  }
#endif
#ifdef PELE_USE_RADIATION
  if (do_rad_solve) {
    BL_PROFILE_VAR("PeleLM::advance::rad", PLM_RAD);
    computeRadSource(a_timestamp_old);
    BL_PROFILE_VAR_STOP(PLM_RAD);
  }
#endif

  addScalarVarianceSources(a_timestamp_old);

  // User defined external sources
  if (m_user_defined_ext_sources) {
    for (int lev = 0; lev <= finest_level; lev++) {
      auto* ldata_p_old = getLevelDataPtr(lev, a_timestamp_old);
      auto* ldata_p_new = getLevelDataPtr(lev, a_timestamp_new);
      auto& ext_src = m_extSource[lev];
      ProblemSpecificFunctions::modify_ext_sources(
        getTime(lev, a_timestamp_old), m_dt, ldata_p_old->state,
        ldata_p_new->state, ext_src, geom[lev].data(), prob_parm_d);
    }
  }
}
