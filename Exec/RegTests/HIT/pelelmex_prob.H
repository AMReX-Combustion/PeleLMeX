#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>
#include <AMReX_REAL.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMFData.H>
#include <PelePhysics.H>
#include <Utilities.H>

using namespace amrex::literals;

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

struct MyProbParm : public ProbParmDefault
{
  // for forcing only
  bool doing_turbulent_forcing = false;
  amrex::Real urms = 1.0;

  // for decay only
  amrex::Real P_mean = 101325.0;
  amrex::Real T_mean = 300.0;

  int input_resolution = 0;
  amrex::Real urms0 = 1.0;
  amrex::Real uin_norm = 1.0;

  amrex::Real Linput = 0.0;

  amrex::Real* d_uinput = nullptr;
  amrex::Real* d_vinput = nullptr;
  amrex::Real* d_winput = nullptr;
  amrex::Real* d_xarray = nullptr;
  amrex::Real* d_xdiff = nullptr;
};

struct MyProblemSpecificFunctions : public DefaultProblemSpecificFunctions
{
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void initdata(
    int i,
    int j,
    int k,
    int is_incompressible,
    amrex::Array4<amrex::Real> const& state,
    amrex::Array4<amrex::Real> const& /*aux*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    constexpr amrex::Real TwoPi = 3.14159265358979323846264338327950288 * 2.0;

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();

    amrex::Real x[3] = {
      prob_lo[0] + static_cast<amrex::Real>(i + 0.5_rt) * dx[0] * 100.0_rt,
      prob_lo[1] + static_cast<amrex::Real>(j + 0.5_rt) * dx[1] * 100.0_rt,
      prob_lo[2] + static_cast<amrex::Real>(k + 0.5_rt) * dx[2] * 100.0_rt};

    if (prob_parm.doing_turbulent_forcing) {
      AMREX_ALWAYS_ASSERT(is_incompressible != 0);
      AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == 3);

      // Fill Velocity
      const amrex::Real tenth_urms = 0.1 * prob_parm.urms;
      const amrex::Real x_norm =
        (x[0] - prob_lo[0]) / (prob_hi[0] - prob_lo[0]) / 100.0;
      const amrex::Real y_norm =
        (x[1] - prob_lo[1]) / (prob_hi[1] - prob_lo[1]) / 100.0;
      const amrex::Real z_norm =
        (x[2] - prob_lo[2]) / (prob_hi[2] - prob_lo[2]) / 100.0;
      state(i, j, k, VELX) = tenth_urms * std::cos(3 * TwoPi * y_norm) *
                             std::cos(3 * TwoPi * z_norm);
      state(i, j, k, VELY) = tenth_urms * std::cos(3 * TwoPi * x_norm) *
                             std::cos(3 * TwoPi * z_norm);
      state(i, j, k, VELZ) = tenth_urms * std::cos(3 * TwoPi * x_norm) *
                             std::cos(3 * TwoPi * y_norm);
    } else { // doing HIT Decay
      // Fill in the velocities
      amrex::Real uinterp[3] = {0.0};

      // Interpolation factors
      amrex::Real mod[3] = {0.0};
      int idx[3] = {0};
      int idxp1[3] = {0};
      amrex::Real slp[3] = {0.0};
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        mod[idim] = std::fmod(x[idim], prob_parm.Linput);
        locate(
          prob_parm.d_xarray, prob_parm.input_resolution, mod[idim], idx[idim]);
        idxp1[idim] = (idx[idim] + 1) % prob_parm.input_resolution;
        slp[idim] = (mod[idim] - prob_parm.d_xarray[idx[idim]]) /
                    prob_parm.d_xdiff[idx[idim]];
      }

      int inSize = prob_parm.input_resolution;

      const amrex::Real f0 = (1 - slp[0]) * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f1 = slp[0] * (1 - slp[1]) * (1 - slp[2]);
      const amrex::Real f2 = (1 - slp[0]) * slp[1] * (1 - slp[2]);
      const amrex::Real f3 = (1 - slp[0]) * (1 - slp[1]) * slp[2];
      const amrex::Real f4 = slp[0] * (1 - slp[1]) * slp[2];
      const amrex::Real f5 = (1 - slp[0]) * slp[1] * slp[2];
      const amrex::Real f6 = slp[0] * slp[1] * (1 - slp[2]);
      const amrex::Real f7 = slp[0] * slp[1] * slp[2];

      uinterp[0] =
        prob_parm.d_uinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 +
        prob_parm.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] *
          f1 +
        prob_parm.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f2 +
        prob_parm.d_uinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f3 +
        prob_parm.d_uinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f4 +
        prob_parm.d_uinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f5 +
        prob_parm.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f6 +
        prob_parm.d_uinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f7;

      uinterp[1] =
        prob_parm.d_vinput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 +
        prob_parm.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] *
          f1 +
        prob_parm.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f2 +
        prob_parm.d_vinput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f3 +
        prob_parm.d_vinput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f4 +
        prob_parm.d_vinput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f5 +
        prob_parm.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f6 +
        prob_parm.d_vinput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f7;

      uinterp[2] =
        prob_parm.d_winput[idx[0] + inSize * (idx[1] + inSize * idx[2])] * f0 +
        prob_parm.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idx[2])] *
          f1 +
        prob_parm.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f2 +
        prob_parm.d_winput[idx[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f3 +
        prob_parm.d_winput[idxp1[0] + inSize * (idx[1] + inSize * idxp1[2])] *
          f4 +
        prob_parm.d_winput[idx[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f5 +
        prob_parm.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idx[2])] *
          f6 +
        prob_parm.d_winput[idxp1[0] + inSize * (idxp1[1] + inSize * idxp1[2])] *
          f7;

      //
      // Fill Velocity
      //
      state(i, j, k, VELX) = uinterp[0];
      state(i, j, k, VELY) = uinterp[1];
      state(i, j, k, VELZ) = uinterp[2];

      if (is_incompressible != 0) {
        return;
      }

      state(i, j, k, TEMP) = prob_parm.T_mean;

      amrex::Real Yt[NUM_SPECIES];
      Yt[O2_ID] = 0.233;
      Yt[N2_ID] = 0.767;

      amrex::Real rho_cgs, P_cgs;
      P_cgs = m2c::P(prob_parm.P_mean);

      auto eos = pele::physics::PhysicsType::eos();
      eos.PYT2R(P_cgs, Yt, state(i, j, k, TEMP), rho_cgs);
      state(i, j, k, DENSITY) = c2m::Rho(rho_cgs);

      eos.TY2H(state(i, j, k, TEMP), Yt, state(i, j, k, RHOH));
      state(i, j, k, RHOH) =
        state(i, j, k, DENSITY) * c2m::H(state(i, j, k, RHOH));

      for (int n = 0; n < NUM_SPECIES; n++) {
        state(i, j, k, FIRSTSPEC + n) = Yt[n] * state(i, j, k, DENSITY);
      }
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
