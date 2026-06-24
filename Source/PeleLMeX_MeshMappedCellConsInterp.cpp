#include <PeleLMeX_MeshMappedCellConsInterp.H>

#include <AMReX_MultiFab.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include <cmath>

namespace {
// Defensive guard: every per-cell multiplication checks that all
// inputs are finite and that fac/detJ are not near-zero, so an
// isolated bad cell can't trigger FE_INVALID under
// amrex.fpe_trap_invalid = 1.  When the guard fails, the cell falls
// back to identity (passes velocity unchanged).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool
guarded_inputs_ok(amrex::Real v, amrex::Real m, amrex::Real d) noexcept
{
  constexpr amrex::Real eps = amrex::Real(1.0e-30);
  return std::isfinite(v) && std::isfinite(m) && std::isfinite(d) &&
         std::abs(m) > eps && std::abs(d) > eps;
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
  AMREX_ASSERT_WITH_MESSAGE(
    ncomp == AMREX_SPACEDIM, "MeshMappedCellConsInterp supports only "
                             "AMREX_SPACEDIM-component velocity data");
  AMREX_ASSERT(m_detJ_crse && m_fac_crse && m_detJ_fine && m_fac_fine);

  // FillPatchTwoLevels passes crsemf as an AMReX-internal temp whose
  // BoxArray is the fine BA coarsened (with the fine DM), so ParallelCopy
  // detJ/fac onto crsemf's layout; setVal(1.0) seeds identity for any
  // cell ParallelCopy leaves untouched.
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

  // Pre-weight: u_xi_i = v_phys_i * detJ / fac_i  (cell-wise, per direction).
  // Standard mf_cell_cons_interp on u_xi IS volume-conservative since Xi-space
  // cells are uniform.
  amrex::MultiFab crse_w(
    crsemf.boxArray(), crsemf.DistributionMap(), ncomp, /*nghost=*/0,
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
        const amrex::Real v = v_arr[box_no](i, j, k, ccomp_l + n);
        const amrex::Real dJ = detJ_arr[box_no](i, j, k);
        const amrex::Real fac = fac_arr[box_no](i, j, k, n);
        if (guarded_inputs_ok(v, fac, dJ)) {
          w_arr[box_no](i, j, k, n) = v * dJ / fac;
        } else {
          w_arr[box_no](i, j, k, n) = std::isfinite(v) ? v : amrex::Real(0.0);
        }
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

  // Unweight: v_phys_i = u_xi_i * fac_i / detJ.  Iterate finemf directly
  // (rather than crse_w + refining indices) so the loop body only touches
  // cells inside finemf[box_no].box() -- critical in 3D where the
  // refinement of crse_w's box does not always match finemf's box, and
  // out-of-bounds Array4 indexing manifests as a bus error.
  {
    auto const& detJ_arr = detJ_fine_local.const_arrays();
    auto const& fac_arr = fac_fine_local.const_arrays();
    auto const& v_arr = finemf.arrays();
    const int fcomp_l = fcomp;
    const amrex::Box dest_dom = dest_domain;
    amrex::ParallelFor(
      finemf, amrex::IntVect(0), ncomp,
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k, int n) noexcept {
        if (dest_dom.contains(amrex::IntVect(AMREX_D_DECL(i, j, k)))) {
          const amrex::Real v = v_arr[box_no](i, j, k, fcomp_l + n);
          const amrex::Real dJ = detJ_arr[box_no](i, j, k);
          const amrex::Real fac = fac_arr[box_no](i, j, k, n);
          if (guarded_inputs_ok(v, dJ, fac)) {
            v_arr[box_no](i, j, k, fcomp_l + n) = v * fac / dJ;
          } else if (!std::isfinite(v)) {
            v_arr[box_no](i, j, k, fcomp_l + n) = amrex::Real(0.0);
          }
          // else: leave v unchanged (identity fallback)
        }
      });
    amrex::Gpu::streamSynchronize();
  }
}
