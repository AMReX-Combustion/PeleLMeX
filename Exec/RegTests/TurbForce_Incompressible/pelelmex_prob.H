#ifndef PELELM_PROB_H
#define PELELM_PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_SPACE.H>

#include <PeleLMeX_ProblemSpecificFunctions.H>
#include <PeleLMeX_Index.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

struct MyProbParm : public ProbParmDefault
{
  amrex::Real urms = 1.0;
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
    AMREX_ALWAYS_ASSERT(is_incompressible != 0);
    AMREX_ALWAYS_ASSERT(AMREX_SPACEDIM == 3);
    constexpr amrex::Real TwoPi = 3.14159265358979323846264338327950288 * 2.0;

    const amrex::Real* prob_lo = geomdata.ProbLo();
    const amrex::Real* prob_hi = geomdata.ProbHi();
    const amrex::Real* dx = geomdata.CellSize();

    AMREX_D_TERM(
      const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
      , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
      , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

    // Fill Velocity
    const amrex::Real tenth_urms = 0.1 * prob_parm.urms;
    const amrex::Real x_norm = (x - prob_lo[0]) / (prob_hi[0] - prob_lo[0]);
    const amrex::Real y_norm = (y - prob_lo[1]) / (prob_hi[1] - prob_lo[1]);
    const amrex::Real z_norm = (z - prob_lo[2]) / (prob_hi[2] - prob_lo[2]);
    state(i, j, k, VELX) =
      tenth_urms * std::cos(3 * TwoPi * y_norm) * std::cos(3 * TwoPi * z_norm);
    state(i, j, k, VELY) =
      tenth_urms * std::cos(3 * TwoPi * x_norm) * std::cos(3 * TwoPi * z_norm);
    state(i, j, k, VELZ) =
      tenth_urms * std::cos(3 * TwoPi * x_norm) * std::cos(3 * TwoPi * y_norm);
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real* /*s_ext[NVAR]*/,
    const int /*idira*/,
    const int /*sgna*/,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdataa*/,
    MyProbParm const& /*prob_parma*/,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;
#endif
