#ifndef PELELM_BCFILL
#define PELELM_BCFILL
#include <PeleLMeX_Index.H>
#include <pelelmex_prob.H>
#include <PMFData.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Geometry.H>

using namespace amrex;

template <int start_fill, int n_fill>
struct PeleLMCCFillExtDirState
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;
  int m_do_turbInflow = 0;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirState(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux,
    int do_turbInflow = 0)
    : lprobparm(a_prob_parm),
      lpmfdata(a_pmf_data),
      m_nAux(a_nAux),
      m_do_turbInflow(do_turbInflow)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& state,
    const int dcomp,
    const int numcomp,
    amrex::GeometryData const& geom,
    const amrex::Real time,
    const amrex::BCRec* bcr,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {

    // Get geometry data
    const int* domlo = geom.Domain().loVect();
    const int* domhi = geom.Domain().hiVect();
    const amrex::Real* prob_lo = geom.ProbLo();
    const amrex::Real* dx = geom.CellSize();
    const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
      prob_lo[0] + (iv[0] + 0.5) * dx[0], prob_lo[1] + (iv[1] + 0.5) * dx[1],
      prob_lo[2] + (iv[2] + 0.5) * dx[2])};

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

      // Low
      for (int n = 0; n < std::min(numcomp, n_fill); n++) {

        // bcnormal handles all the state components at once
        amrex::Real s_ext[NVAR] = {0.0};

        // Get current state component BC
        const int* bc = bcr[n].data();

        if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

          // If using TurbInflow, pass in the turbulent data. User can overwrite
          // if needed
          if (m_do_turbInflow != 0) {
            if (n >= VELX && n < DENSITY) {
              if (!isnan(
                    state(iv, VELX + n))) { // During fillcoarsepatch, get here
                                            // with nan. TODO: find a better fix
                s_ext[VELX + n] = state(iv, VELX + n);
              }
            }
          }
          // bcnormal() is defined in pelelmex_prob.H in problem directory in
          // /Exec
          bcnormal(x, m_nAux, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);
          state(iv, n) = s_ext[n];
        }
      }

      // High
      for (int n = 0; n < std::min(numcomp, n_fill); n++) {

        // bcnormal handles all the state components at once
        amrex::Real s_ext[NVAR] = {0.0};

        // Get current state component BC
        const int* bc = bcr[n].data();

        if (
          (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
          (iv[idir] > domhi[idir])) {

          // If using TurbInflow, pass in the turbulent data. User can overwrite
          // if needed
          if (m_do_turbInflow != 0) {
            if (start_fill + n >= VELX && start_fill + n < DENSITY) {
              if (!isnan(
                    state(iv, VELX + n))) { // During fillcoarsepatch, get here
                                            // with nan. TODO: find a better fix
                s_ext[VELX + n] = state(iv, VELX + n);
              }
            }
          }
          // bcnormal() is defined in pelelmex_prob.H in problem directory in
          // /Exec
          bcnormal(
            x, m_nAux, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);
          state(iv, dcomp + n) = s_ext[start_fill + n];
        }
      }
    }
  }
};

struct PeleLMCCFillExtDirDummy
{

  ProbParm const* lprobparm;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirDummy(ProbParm const* a_prob_parm, int a_nAux)
    : lprobparm(a_prob_parm), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dummy*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
  }
};

//
// A dummy function because FillPatch requires something to exist for filling
// dirichlet boundary conditions, even if we know we cannot have an ext_dir BC.
// u_mac BCs are only either periodic (amrex::BCType::int_dir) or first order
// extrapolation (amrex::BCType::foextrap).
//
struct umacFill
{
  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dummy*/,
    const int /*dcomp*/,
    const int numcomp,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* bcr,
    const int bcomp,
    const int /*orig_comp*/) const
  {
    // Abort if this function is expected to fill an ext_dir BC.
    for (int n = bcomp; n < bcomp + numcomp; ++n) {
      const amrex::BCRec& bc = bcr[n];
      if (AMREX_D_TERM(
            bc.lo(0) == amrex::BCType::ext_dir ||
              bc.hi(0) == amrex::BCType::ext_dir,
            || bc.lo(1) == amrex::BCType::ext_dir ||
              bc.hi(1) == amrex::BCType::ext_dir,
            || bc.lo(2) == amrex::BCType::ext_dir ||
              bc.hi(2) == amrex::BCType::ext_dir)) {
        amrex::Abort("PeleLM::umacFill: umac should not have BCType::ext_dir");
      }
    }
  }
};
#endif
