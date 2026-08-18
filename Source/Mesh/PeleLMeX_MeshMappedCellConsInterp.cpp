#include <PeleLMeX_MeshMappedCellConsInterp.H>

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>

namespace {
constexpr amrex::Real meshmap_interp_eps = amrex::Real(1.0e-30);

// Per-component Xi-space weight: velocity -> detJ/fac_i, conserved
// density -> detJ, else -> 1.  Falls back to 1 if the metrics are
// non-finite or near-zero (avoids FE_INVALID under fpe_trap_invalid).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
metric_weight(
  bool is_vel, bool is_dens, amrex::Real detJ, amrex::Real fac) noexcept
{
  if (is_vel) {
    if (
      std::isfinite(detJ) && std::isfinite(fac) &&
      std::abs(detJ) > meshmap_interp_eps &&
      std::abs(fac) > meshmap_interp_eps) {
      return detJ / fac;
    }
    return amrex::Real(1.0);
  }
  if (is_dens) {
    if (std::isfinite(detJ) && std::abs(detJ) > meshmap_interp_eps) {
      return detJ;
    }
    return amrex::Real(1.0);
  }
  return amrex::Real(1.0);
}
} // namespace

void
MeshMappedCellConsInterp::interp(
  amrex::MultiFab const& crsemf,
  int ccomp,
  amrex::MultiFab& finemf,
  int fcomp,
  int ncomp,
  amrex::IntVect const& ng,
  amrex::Geometry const& cgeom,
  amrex::Geometry const& fgeom,
  amrex::Box const& dest_domain,
  amrex::IntVect const& ratio,
  amrex::Vector<amrex::BCRec> const& bcs,
  int bcscomp)
{
  AMREX_ASSERT(m_detJ_crse && m_fac_crse && m_detJ_fine && m_fac_fine);

  // Absolute component = ccomp/fcomp + n (PeleLM fills the whole state
  // from component 0); used to pick the per-component weight.
  const int vel_start = m_vel_start;
  const int vel_ncomp = m_vel_ncomp;
  const int dens_start = m_dens_start;
  const int dens_end = m_dens_end;

  // crsemf is an AMReX temp on the coarsened-fine layout, so ParallelCopy
  // the metrics onto it; setVal(1.0) leaves uncovered cells at identity.
  amrex::MultiFab detJ_crse_local(
    crsemf.boxArray(), crsemf.DistributionMap(), 1, 0, amrex::MFInfo(),
    crsemf.Factory());
  detJ_crse_local.setVal(1.0);
  detJ_crse_local.ParallelCopy(
    *m_detJ_crse, 0, 0, 1, m_detJ_crse->nGrowVect(), amrex::IntVect{0});

  amrex::MultiFab fac_crse_local(
    crsemf.boxArray(), crsemf.DistributionMap(), AMREX_SPACEDIM, 0,
    amrex::MFInfo(), crsemf.Factory());
  fac_crse_local.setVal(1.0);
  fac_crse_local.ParallelCopy(
    *m_fac_crse, 0, 0, AMREX_SPACEDIM, m_fac_crse->nGrowVect(),
    amrex::IntVect{0});

  // Pre-weight the coarse data into Xi-space (v_xi = v * w).  Standard
  // mf_cell_cons_interp on v_xi IS volume-conservative since Xi-space
  // cells are uniform.
  amrex::MultiFab crse_w(
    crsemf.boxArray(), crsemf.DistributionMap(), ncomp, /*ngrow=*/0,
    amrex::MFInfo(), crsemf.Factory());
  {
    auto const& v_arr = crsemf.const_arrays();
    auto const& detJ_arr = detJ_crse_local.const_arrays();
    auto const& fac_arr = fac_crse_local.const_arrays();
    auto const& w_arr = crse_w.arrays();
    const int ccomp_l = ccomp;
    amrex::ParallelFor(
      crse_w, amrex::IntVect(0), ncomp,
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
        const int c = ccomp_l + n;
        const bool is_vel = (c >= vel_start && c < vel_start + vel_ncomp);
        const bool is_dens = (c >= dens_start && c <= dens_end);
        const amrex::Real v = v_arr[box_no](i, j, k, ccomp_l + n);
        const amrex::Real dJ = detJ_arr[box_no](i, j, k);
        const amrex::Real fac =
          is_vel ? fac_arr[box_no](i, j, k, c - vel_start) : amrex::Real(1.0);
        const amrex::Real w = metric_weight(is_vel, is_dens, dJ, fac);
        w_arr[box_no](i, j, k, n) = std::isfinite(v) ? v * w : amrex::Real(0.0);
      });
    amrex::Gpu::streamSynchronize();
  }

  // Standard cell-cons interp on the Xi-space coarse data.
  amrex::Vector<amrex::BCRec> bcs_local(
    bcs.begin() + bcscomp, bcs.begin() + bcscomp + ncomp);
  amrex::mf_cell_cons_interp.interp(
    crse_w, 0, finemf, fcomp, ncomp, ng, cgeom, fgeom, dest_domain, ratio,
    bcs_local, 0);

  // Materialize fine-side metrics on finemf's BoxArray (same pattern,
  // same identity-default rationale).
  amrex::MultiFab detJ_fine_local(
    finemf.boxArray(), finemf.DistributionMap(), 1, ng, amrex::MFInfo(),
    finemf.Factory());
  detJ_fine_local.setVal(1.0);
  detJ_fine_local.ParallelCopy(
    *m_detJ_fine, 0, 0, 1, m_detJ_fine->nGrowVect(), ng);

  amrex::MultiFab fac_fine_local(
    finemf.boxArray(), finemf.DistributionMap(), AMREX_SPACEDIM, ng,
    amrex::MFInfo(), finemf.Factory());
  fac_fine_local.setVal(1.0);
  fac_fine_local.ParallelCopy(
    *m_fac_fine, 0, 0, AMREX_SPACEDIM, m_fac_fine->nGrowVect(), ng);

  // Un-weight back to physical space (v = v_xi / w) over the SAME region
  // the cell-cons interp wrote: finemf valid + ng ghost, clipped to
  // dest_domain.  The ng (coarse/fine) ghost cells carry Xi-space weights
  // under a non-trivial map and must be un-weighted too.  Iterate finemf
  // directly (always in-bounds) -- refining crse_w's box can stray out of
  // bounds in 3D (bus error).
  {
    auto const& detJ_arr = detJ_fine_local.const_arrays();
    auto const& fac_arr = fac_fine_local.const_arrays();
    auto const& v_arr = finemf.arrays();
    const int fcomp_l = fcomp;
    const amrex::Box dest_dom = dest_domain;
    amrex::ParallelFor(
      finemf, ng, ncomp,
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
        if (dest_dom.contains(amrex::IntVect(AMREX_D_DECL(i, j, k)))) {
          const int c = fcomp_l + n;
          const bool is_vel = (c >= vel_start && c < vel_start + vel_ncomp);
          const bool is_dens = (c >= dens_start && c <= dens_end);
          const amrex::Real v = v_arr[box_no](i, j, k, fcomp_l + n);
          const amrex::Real dJ = detJ_arr[box_no](i, j, k);
          const amrex::Real fac =
            is_vel ? fac_arr[box_no](i, j, k, c - vel_start) : amrex::Real(1.0);
          const amrex::Real w = metric_weight(is_vel, is_dens, dJ, fac);
          if (std::isfinite(v)) {
            v_arr[box_no](i, j, k, fcomp_l + n) = v / w;
          } else {
            v_arr[box_no](i, j, k, fcomp_l + n) = amrex::Real(0.0);
          }
        }
      });
    amrex::Gpu::streamSynchronize();
  }
}
