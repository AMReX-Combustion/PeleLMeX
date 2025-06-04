#ifndef PELELM_BCFILL
#define PELELM_BCFILL
#include <PeleLMeX_Index.H>
#include <pelelmex_prob.H>
#include <PMFData.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Geometry.H>

using namespace amrex;

template <typename T>
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
    int do_turbInflow)
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
    const int /*dcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

      // Low
      for (int n = 0; n < std::min(numcomp, NVAR); n++) {

        // bcnormal handles all the state components at once
        amrex::Real s_ext[NVAR] = {0.0};
        // Get current state component BC
        const int* bc = bcr[n].data();

        if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

          // Get the state in the first cell in the domain
          auto iv_in = iv_3d;
          iv_in[idir] = iv_in[idir] + 1;
          const auto s_in = state.cellData(iv_in[0], iv_in[1], iv_in[2]);

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
          T::bcnormal(
            x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);
          state(iv, n) = s_ext[n];
        }
      }

      // High
      for (int n = 0; n < std::min(numcomp, NVAR); n++) {

        // bcnormal handles all the state components at once
        amrex::Real s_ext[NVAR] = {0.0};

        // Get current state component BC
        const int* bc = bcr[n].data();

        if (
          (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
          (iv[idir] > domhi[idir])) {

          // Get the state in the first cell in the domain
          auto iv_in = iv_3d;
          iv_in[idir] = iv_in[idir] - 1;
          const auto s_in = state.cellData(iv_in[0], iv_in[1], iv_in[2]);

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
          T::bcnormal(
            x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);
          state(iv, n) = s_ext[n];
        }
      }
    }
  }
};

template <typename T>
struct PeleLMCCFillExtDirAux
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirAux(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& aux,
    const int /*dcomp*/,
    const int /*numcomp*/,
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

#if AMREX_SPACEDIM == 2
    auto aux_ext = aux.cellData(iv[0], iv[1], 0);
#else
    auto aux_ext = aux.cellData(iv[0], iv[1], iv[2]);
#endif

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {

      // Low
      for (int n = 0; n < m_nAux; n++) {

        // Get current state component BC
        const int* bc = bcr[n].data();

        if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {
          // bcnormal_aux() is defined in pelelmex_prob.H in problem directory
          // in /Exec
          T::bcnormal_aux(
            x, aux_ext, idir, 1, time, geom, *lprobparm, lpmfdata);
        }
      }
      // High
      for (int n = 0; n < m_nAux; n++) {

        // Get current state component BC
        const int* bc = bcr[n].data();

        if (
          (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
          (iv[idir] > domhi[idir])) {
          T::bcnormal_aux(
            x, aux_ext, idir, -1, time, geom, *lprobparm, lpmfdata);
        }
      }
    }
  }
};

template <typename T>
struct PeleLMCCFillExtDirSpec
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirSpec(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& rhoY,
    const int dcomp,
    const int /*numcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    for (int n = 0; n < NUM_SPECIES; n++) {

      // Get the species component BC
      const int* bc = bcr[n].data();

      // bcnormal handles all the state components at once
      amrex::Real s_ext[NVAR] = {0.0};

      for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
        // Low
        if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

          // Get the state on the first cell in the domain
          auto iv_in = iv_3d;
          iv_in[idir] = iv_in[idir] + 1;
          const auto s_in = rhoY.cellData(iv_in[0], iv_in[1], iv_in[2]);

          // bcnormal() is defined in pelelmex_prob.H in problem directory in
          // /Exec
          T::bcnormal(
            x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

          rhoY(iv, dcomp + n) = s_ext[FIRSTSPEC + n];

          // High
        } else if (
          (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
          (iv[idir] > domhi[idir])) {

          // Get the state on the first cell in the domain
          auto iv_in = iv_3d;
          iv_in[idir] = iv_in[idir] - 1;
          const auto s_in = rhoY.cellData(iv_in[0], iv_in[1], iv_in[2]);

          // bcnormal() is defined in pelelmex_prob.H in problem directory in
          // /Exec
          T::bcnormal(
            x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

          rhoY(iv, dcomp + n) = s_ext[FIRSTSPEC + n];
        }
      }
    }
  }
};

template <typename T>
struct PeleLMCCFillExtDirDens
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirDens(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& state,
    const int dcomp,
    const int /*numcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    // Get the first species component BC
    const int* bc = bcr->data();

    // bcnormal handles all the state components at once
    amrex::Real s_ext[NVAR] = {0.0};

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      // Low
      if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] + 1;
        const auto s_in = state.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

        state(iv, dcomp) = s_ext[DENSITY];

        // High
      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
        (iv[idir] > domhi[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] - 1;
        const auto s_in = state.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

        state(iv, dcomp) = s_ext[DENSITY];
      }
    }
  }
};

template <typename T>
struct PeleLMCCFillExtDirRhoH
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirRhoH(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& rhoH,
    const int /*dcomp*/,
    const int /*numcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    // Get the first species component BC
    const int* bc = bcr->data();

    // bcnormal handles all the state components at once
    amrex::Real s_ext[NVAR] = {0.0};

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      // Low
      if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] + 1;
        const auto s_in = rhoH.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

        rhoH(iv) = s_ext[RHOH];

        // High
      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
        (iv[idir] > domhi[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] - 1;
        const auto s_in = rhoH.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

        rhoH(iv) = s_ext[RHOH];
      }
    }
  }
};

template <typename T>
struct PeleLMCCFillExtDirTemp
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirTemp(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& state,
    const int dcomp,
    const int /*numcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    // Get the first species component BC
    const int* bc = bcr->data();

    // bcnormal handles all the state components at once
    amrex::Real s_ext[NVAR] = {0.0};
    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      // Low
      if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] + 1;
        const auto s_in = state.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

        state(iv, dcomp) = s_ext[TEMP];

        // High
      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
        (iv[idir] > domhi[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] - 1;
        const auto s_in = state.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

        state(iv, dcomp) = s_ext[TEMP];
      }
    }
  }
};

#ifdef PELE_USE_PLASMA
template <typename T>
struct PeleLMCCFillExtDirnE
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirnE(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& nE,
    const int /*dcomp*/,
    const int /*numcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    // Get the first species component BC
    const int* bc = bcr->data();

    // bcnormal handles all the state components at once
    amrex::Real s_ext[NVAR] = {0.0};

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      // Low
      if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] + 1;
        const auto s_in = nE.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

        nE(iv) = s_ext[NE];

        // High
      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
        (iv[idir] > domhi[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] - 1;
        const auto s_in = nE.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

        nE(iv) = s_ext[NE];
      }
    }
  }
};
#endif

#ifdef PELE_USE_PLASMA
template <typename T>
struct PeleLMCCFillExtDirPhiV
{

  ProbParm const* lprobparm;
  pele::physics::PMF::PmfData::DataContainer const* lpmfdata;
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirPhiV(
    ProbParm const* a_prob_parm,
    pele::physics::PMF::PmfData::DataContainer const* a_pmf_data,
    int a_nAux)
    : lprobparm(a_prob_parm), lpmfdata(a_pmf_data), m_nAux(a_nAux)
  {
  }

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& iv,
    amrex::Array4<amrex::Real> const& phiV,
    const int dcomp,
    const int /*numcomp*/,
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

    // Get a 3D version of iv to access first cell in domain
    auto iv_3d = amrex::IntVectND<3>(0);
    AMREX_D_TERM(iv_3d[0] = iv[0];, iv_3d[1] = iv[1];, iv_3d[2] = iv[2];);

    // Get the first species component BC
    const int* bc = bcr->data();

    // bcnormal handles all the state components at once
    amrex::Real s_ext[NVAR] = {0.0};

    for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
      // Low
      if ((bc[idir] == amrex::BCType::ext_dir) and (iv[idir] < domlo[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] + 1;
        const auto s_in = phiV.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, 1, time, geom, *lprobparm, lpmfdata);

        phiV(iv, dcomp) = s_ext[PHIV];

        // High
      } else if (
        (bc[idir + AMREX_SPACEDIM] == amrex::BCType::ext_dir) and
        (iv[idir] > domhi[idir])) {

        // Get the state on the first cell in the domain
        auto iv_in = iv_3d;
        iv_in[idir] = iv_in[idir] - 1;
        const auto s_in = phiV.cellData(iv_in[0], iv_in[1], iv_in[2]);

        // bcnormal() is defined in pelelmex_prob.H in problem directory in
        // /Exec
        T::bcnormal(x, s_in, s_ext, idir, -1, time, geom, *lprobparm, lpmfdata);

        phiV(iv, dcomp) = s_ext[PHIV];
      }
    }
  }
};
#endif

struct PeleLMCCFillExtDirDummy
{
  const int m_nAux;

  AMREX_GPU_HOST
  constexpr PeleLMCCFillExtDirDummy(int a_nAux) : m_nAux(a_nAux) {}

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
