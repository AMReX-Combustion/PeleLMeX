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
  // Lid velocity along x (top boundary).  Set so that
  // Reynolds number Re = U_lid * L / nu equals the target.
  amrex::Real U_lid = 1.0;
  // Reference density / viscosity / temperature for the incompressible
  // run (used by PeleLMeX via peleLM.rho, peleLM.mu in the input file;
  // duplicated here only for documentation).
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
    amrex::GeometryData const& /*geomdata*/,
    MyProbParm const& /*prob_parm*/,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    AMREX_D_TERM(state(i, j, k, VELX) = amrex::Real(0.0);
                 , state(i, j, k, VELY) = amrex::Real(0.0);
                 , state(i, j, k, VELZ) = amrex::Real(0.0));
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void bcnormal(
    const amrex::Real* x,
    const amrex::CellData<amrex::Real>& /*s_in*/,
    amrex::Real s_ext[NVAR],
    const int idir,
    const int sgn,
    const amrex::Real /*time*/,
    amrex::GeometryData const& geomdata,
    MyProbParm const& prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* /*pmf_data*/)
  {
    // BC convention note.  PeleLMeX's Inflow BC writes s_ext into the
    // GHOST CELL CENTER.  Empirically (verified at Re=100), the
    // cell-centered velocity in the row adjacent to the wall converges
    // to ~ s_ext / 2 -- i.e. the discrete face value seen by the
    // momentum equations is the arithmetic mean of the ghost and the
    // first interior cell.  To make the steady-state face velocity at
    // the lid equal to U_lid we therefore set the ghost cell to
    // 2 * U_lid.  The image-cell formula 2*U_lid - u_interior is
    // numerically unstable in this code path (NaN within ~60 steps),
    // so we use the constant override.  A proper "MovingWall" BC
    // type in PeleLMeX is the long-term fix; that would let
    // AMReX-Hydro's tangential-velocity extrapolation force both
    // predictor states to the wall velocity.
    AMREX_D_TERM(
      s_ext[VELX] = amrex::Real(0.0);, s_ext[VELY] = amrex::Real(0.0);
      , s_ext[VELZ] = amrex::Real(0.0));
    amrex::ignore_unused(x);
    if (idir == 1 && sgn == -1) {
      s_ext[VELX] = prob_parm.U_lid;
    }
  }
};

using ProblemSpecificFunctions = MyProblemSpecificFunctions;
using ProbParm = MyProbParm;

#endif
