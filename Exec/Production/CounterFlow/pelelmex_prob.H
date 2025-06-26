#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMFData.H>
#include <PelePhysics.H>

struct MyProbParm : public ProbParmDefault
{
  // amrex::Real T_ox     = 298.0;
  // amrex::Real T_fuel   = 298.0;
  // amrex::Real T_inert   = 298.0;
  amrex::Real T_ox = 450.0;
  amrex::Real T_fuel = 450.0;
  amrex::Real T_inert = 450.0;
  amrex::Real Y_O2_ox = 0.233;
  amrex::Real Y_N2_ox = 0.767;
  amrex::Real Y_N2_fuel = 0.0;
  amrex::Real Y_fuel = 0.0;
  amrex::Real Y_H2O_ign = 0.0;
  amrex::Real Y_CO2_ign = 0.0;
  amrex::Real mdot_ox = 0.01;
  amrex::Real mdot_fuel = 0.01;
  amrex::Real pertmag = 0.0004;
  amrex::Real inertVel_ox = 0.01;
  amrex::Real inertVel_fuel = 0.01;
  amrex::Real jetRadius = 0.004;
  amrex::Real inertRadius = 0.005;

  // Misc
  int do_ignit = 0;
  amrex::Real ignitSphereRad = 0.003;
  amrex::Real ignitT = 1800.0;
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata(
    int i,
    int j,
    int k,
    int /*is_incompressible*/,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();

    AMREX_D_TERM(
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
      , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

    AMREX_D_TERM(
      const amrex::Real Lx = prob_hi[0] - prob_lo[0];
      , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
      , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    AMREX_D_TERM(
      const amrex::Real xc = prob_lo[0] + Lx / 2.0;
      , const amrex::Real yc = prob_lo[1] + Ly / 2.0;
      , const amrex::Real zc = prob_lo[2] + Lz / 2.0);

    constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;

    auto eos = pele::physics::PhysicsType::eos();
    amrex::Real massfrac[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_cold[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_hot[NUM_SPECIES] = {0.0};

    massfrac_cold[NC12H26_ID] = 0.05;
    massfrac_cold[H2O_ID] = 0.0;
    massfrac_cold[CO2_ID] = 0.0;
    massfrac_cold[O2_ID] = 0.233 * (1.0 - massfrac_cold[NC12H26_ID]);
    ;
    massfrac_cold[N2_ID] =
      1.0 - massfrac_cold[NC12H26_ID] - massfrac_cold[O2_ID];

    state(i, j, k, TEMP) = prob_parm.T_inert;

    if (prob_parm.do_ignit) {
      amrex::Real radiusSq = AMREX_D_TERM(
        (x - xc) * (x - xc), +(y - yc) * (y - yc), +(z - zc) * (z - zc));
      amrex::Real radius = std::sqrt(radiusSq);
      amrex::Real mixingWidth = 0.1 * prob_parm.ignitSphereRad;
      amrex::Real mixingFunction =
        0.5 *
        (1.0 + std::tanh((prob_parm.ignitSphereRad - radius) / mixingWidth));
      state(i, j, k, TEMP) = mixingFunction * prob_parm.ignitT +
                             (1.0 - mixingFunction) * prob_parm.T_inert;
      massfrac_hot[NC12H26_ID] = 0.00;
      massfrac_hot[H2O_ID] = prob_parm.Y_H2O_ign;
      massfrac_hot[CO2_ID] = prob_parm.Y_CO2_ign;
      massfrac_hot[O2_ID] =
        0.233 * (1.0 - massfrac_hot[H2O_ID] - massfrac_hot[CO2_ID]);
      massfrac_hot[N2_ID] =
        1.0 - massfrac_hot[H2O_ID] - massfrac_hot[CO2_ID] - massfrac_hot[O2_ID];

      massfrac[NC12H26_ID] = mixingFunction * massfrac_hot[NC12H26_ID] +
                             (1.0 - mixingFunction) * massfrac_cold[NC12H26_ID];
      massfrac[H2O_ID] = mixingFunction * massfrac_hot[H2O_ID] +
                         (1.0 - mixingFunction) * massfrac_cold[H2O_ID];
      massfrac[CO2_ID] = mixingFunction * massfrac_hot[CO2_ID] +
                         (1.0 - mixingFunction) * massfrac_cold[CO2_ID];
      massfrac[O2_ID] = mixingFunction * massfrac_hot[O2_ID] +
                        (1.0 - mixingFunction) * massfrac_cold[O2_ID];
      massfrac[N2_ID] = mixingFunction * massfrac_hot[N2_ID] +
                        (1.0 - mixingFunction) * massfrac_cold[N2_ID];
    } else {
      massfrac[NC12H26_ID] = massfrac_cold[NC12H26_ID];
      massfrac[H2O_ID] = massfrac_cold[H2O_ID];
      massfrac[CO2_ID] = massfrac_cold[CO2_ID];
      massfrac[O2_ID] = massfrac_cold[O2_ID];
      massfrac[N2_ID] = massfrac_cold[N2_ID];
    }

    AMREX_D_TERM(state(i, j, k, VELX) = 0.0;, state(i, j, k, VELY) = 0.0;
                 , state(i, j, k, VELZ) = 0.0);

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;

    // Density
    amrex::Real rho_cgs = 0.0;
    eos.PYT2R(P_cgs, massfrac, state(i, j, k, TEMP), rho_cgs);
    state(i, j, k, DENSITY) = rho_cgs * 1.0e3;

    // Enthalpy
    amrex::Real h_cgs = 0.0;
    eos.TY2H(state(i, j, k, TEMP), massfrac, h_cgs);
    state(i, j, k, RHOH) = h_cgs * 1.0e-4 * state(i, j, k, DENSITY);

    // Species mass
    for (int n = 0; n < NUM_SPECIES; n++) {
      state(i, j, k, FIRSTSPEC + n) = massfrac[n] * state(i, j, k, DENSITY);
    }
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real x[AMREX_SPACEDIM],
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int idir,
    const int sgn,
    const amrex::Real time,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();

    AMREX_D_TERM(
      const amrex::Real Lx = prob_hi[0] - prob_lo[0];
      , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
      , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

    amrex::Real splitx = prob_lo[0] + 0.5 * Lx;
    amrex::Real splity = prob_lo[1] + 0.5 * Ly;

    amrex::Real massfrac_ox[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_fuel[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_inert[NUM_SPECIES] = {0.0};
    amrex::Real massfrac_mix[NUM_SPECIES] = {0.0};

    massfrac_ox[N2_ID] = prob_parm.Y_N2_ox;
    massfrac_ox[O2_ID] = prob_parm.Y_O2_ox;
    massfrac_fuel[NC12H26_ID] = prob_parm.Y_fuel; // Replace ID with gas fuel
    massfrac_fuel[N2_ID] = prob_parm.Y_N2_fuel;
    massfrac_inert[N2_ID] = 1.0;

    auto eos = pele::physics::PhysicsType::eos();

    for (int n = 0; n < NVAR; n++) {
      s_ext[n] = 0.0;
    }

    // Get zone
    int zone = -1;
    if (x[0] < splitx) {
      zone = 1; // Oxidizer side jet
    } else {
      zone = 2; // Fuel side jet
    }
    amrex::Real radius = std::abs(x[1] - splity);
    amrex::Real mixingWidth = 0.1 * prob_parm.jetRadius;
    amrex::Real mixingFunction =
      0.5 * (1.0 + std::tanh((prob_parm.jetRadius - radius) / mixingWidth));

    amrex::Real P_cgs = prob_parm.P_mean * 10.0;
    amrex::Real rho_cgs, RhoH_temp;
    amrex::Real rhoOx = -1.0;
    {
      eos.PYT2R(P_cgs, massfrac_ox, prob_parm.T_ox, rhoOx);
      rhoOx *= 1.0e3;
    }
    amrex::Real rhoFuel = -1.0;
    {
      eos.PYT2R(P_cgs, massfrac_fuel, prob_parm.T_fuel, rhoFuel);
      rhoFuel *= 1.0e3;
    }
    if (zone == 1) {
      s_ext[TEMP] = prob_parm.T_ox * mixingFunction +
                    prob_parm.T_inert * (1.0 - mixingFunction);
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_mix[n] = massfrac_ox[n] * mixingFunction +
                          massfrac_inert[n] * (1.0 - mixingFunction);
      }
      eos.PYT2R(P_cgs, massfrac_mix, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac_mix, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac_mix[n] * s_ext[DENSITY];
      }
    } else if (zone == 2) {
      s_ext[TEMP] = prob_parm.T_fuel * mixingFunction +
                    prob_parm.T_inert * (1.0 - mixingFunction);
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac_mix[n] = massfrac_fuel[n] * mixingFunction +
                          massfrac_inert[n] * (1.0 - mixingFunction);
      }
      eos.PYT2R(P_cgs, massfrac_mix, s_ext[TEMP], rho_cgs);
      s_ext[DENSITY] = rho_cgs * 1.0e3;

      eos.TY2H(s_ext[TEMP], massfrac_mix, RhoH_temp);
      s_ext[RHOH] = RhoH_temp * 1.0e-4 * s_ext[DENSITY];

      for (int n = 0; n < NUM_SPECIES; n++) {
        s_ext[FIRSTSPEC + n] = massfrac_mix[n] * s_ext[DENSITY];
      }
    }

    s_ext[VELY] = 0.0;
    if (zone == 1) {
      amrex::Real vel = (prob_parm.mdot_ox / rhoOx) * mixingFunction +
                        prob_parm.inertVel_ox /* (1.0 + (radius/0.025))*/ *
                          (1.0 - mixingFunction);
      s_ext[VELX] = vel;
    } else if (zone == 2) {
      amrex::Real vel = (prob_parm.mdot_fuel / rhoFuel) * mixingFunction +
                        prob_parm.inertVel_fuel /* (1.0 + (radius/0.025))*/ *
                          (1.0 - mixingFunction);
      s_ext[VELX] = -vel;
    }
  }

  static void patchFlowVariables(
    const amrex::Geometry& geom,
    MyProbParm const& lprobparm,
    amrex::MultiFab& a_mf)
  {

    amrex::Print() << "\nPatching flow variables..";

    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* prob_hi = geom.ProbHi();
    const amrex::Real* dx = geom.CellSize();

    for (amrex::MFIter mfi(a_mf, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.tilebox();

      auto const& rho_arr = a_mf.array(mfi, DENSITY);
      auto const& rhoY_arr = a_mf.array(mfi, FIRSTSPEC);
      auto const& rhoH_arr = a_mf.array(mfi, RHOH);
      auto const& temp_arr = a_mf.array(mfi, TEMP);
      amrex::Real massfrac[NUM_SPECIES] = {0.0};

      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          auto eos = pele::physics::PhysicsType::eos();
          amrex::Real massfrac[NUM_SPECIES] = {0.0};
          amrex::Real sumYs = 0.0;

          for (int n = 0; n < NUM_SPECIES; n++) {
            massfrac[n] = rhoY_arr(i, j, k, n);
#ifdef N2_ID
            if (n != N2_ID) {
              sumYs += massfrac[n];
            }
#endif
          }
#ifdef N2_ID
          massfrac[N2_ID] = 1.0 - sumYs;
#endif

          AMREX_D_TERM(
            const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
            , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
            , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

          AMREX_D_TERM(
            const amrex::Real Lx = prob_hi[0] - prob_lo[0];
            , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
            , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

          AMREX_D_TERM(
            const amrex::Real xc = prob_lo[0] + Lx / 2.0;
            , const amrex::Real yc = prob_lo[1] + Ly / 2.0;
            , const amrex::Real zc = prob_lo[2] + Lz / 2.0);

          amrex::Real radiusSq = AMREX_D_TERM(
            (x - xc) * (x - xc), +(y - yc) * (y - yc), +(z - zc) * (z - zc));

          amrex::Real radius = std::sqrt(radiusSq);
          amrex::Real mixingWidth = 0.1 * lprobparm.ignitSphereRad;
          amrex::Real mixingFunction =
            0.5 * (1.0 + std::tanh(
                           (lprobparm.ignitSphereRad - radius) / mixingWidth));
          temp_arr(i, j, k) = mixingFunction * lprobparm.ignitT +
                              (1.0 - mixingFunction) * lprobparm.T_inert;
        });
    }

    amrex::Print() << "Done\n";
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
