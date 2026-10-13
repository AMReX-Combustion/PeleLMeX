#ifndef PELELMEX_MESHMAP_EVALUATOR_H
#define PELELMEX_MESHMAP_EVALUATOR_H

// The mesh-mapping descriptor lives in PelePhysics
// (Source/Utility/MeshMap/MeshMapEvaluator.H) so that the turbulent-inflow
// reader and the inflow-file generators share the solver's definition of
// the mapping.  PeleLMeX code keeps referring to it by its unqualified name.
#include <MeshMapEvaluator.H>

using MeshMapEvaluator = pele::physics::MeshMapEvaluator;

/**
 * \brief Physical cell-centre coordinate along \p idir for a possibly
 *        out-of-domain index, with periodic directions wrapped.
 *
 * MeshMapEvaluator::x_phys_cc continues the map past the domain end.  For
 * a ghost cell across a periodic boundary the physical position must be
 * the periodic image of the interior cell it represents (x of that cell
 * shifted by the physical period); a non-linear map continued past the
 * end gives a different point, and inflow values written by position at
 * such ghost cells (bcnormal, turbulent inflow) then differ between the
 * two periodic copies of the same cell, which the nodal projection
 * cannot reconcile.  Non-periodic directions keep the continued map,
 * which inflow faces rely on.  Identical to x_phys_cc for indices inside
 * the domain, and bypassed (bit-identical positions) for the identity map,
 * where the continued map already is the periodic image.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
x_phys_cc_periodic(
  MeshMapEvaluator const& mmap,
  int idir,
  int idx,
  amrex::GeometryData const& gd) noexcept
{
  const int lo = gd.Domain().smallEnd(idir);
  const int n = gd.Domain().length(idir);
  if (
    mmap.m_kind != MeshMapEvaluator::Kind::Identity &&
    gd.isPeriodic(idir) != 0 && (idx < lo || idx >= lo + n)) {
    // shift is the number of periods between idx and its image in the domain
    const int k = idx - lo;
    const int shift = (k >= 0) ? k / n : -((-k - 1) / n + 1);
    return mmap.x_phys_cc(idir, idx - shift * n, gd) +
           static_cast<amrex::Real>(shift) * mmap.domain_size_phys(idir, gd);
  }
  return mmap.x_phys_cc(idir, idx, gd);
}

/**
 * \brief Inverse physical spacing for a centred difference at cell \p idx
 *        along \p idir: 2 / (x_{idx+1} - x_{idx-1}), i.e. the factor that
 *        turns a Xi-index central difference into a physical derivative.
 *
 * Without mapping (\p mapped false) returns \p inv_dx_uniform unchanged,
 * so unmapped stencils stay bit-identical to the historical ones; with it,
 * (a_{i+1} - a_{i-1}) * 0.5 * result is the second-order physical derivative
 * on the stretched grid.  The neighbour positions come from
 * x_phys_cc_periodic, so at a periodic seam the spacing is the physical
 * distance to the periodic image (what the metric ghost cells hold), not
 * the map continued past the domain end; in non-periodic directions the
 * continued map is used, consistent with the fill-patched ghost state the
 * stencil reads there.  Positions are analytic evaluations of the map, so
 * any index is admissible.  Used by the gradient derives (vorticity,
 * Q-criterion, enstrophy) and the divU time-step estimate.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
mapped_inv_dx(
  int idir,
  int idx,
  bool mapped,
  amrex::Real inv_dx_uniform,
  MeshMapEvaluator const& mmap,
  amrex::GeometryData const& gd) noexcept
{
  if (!mapped) {
    return inv_dx_uniform;
  }
  return amrex::Real(2.0) / (x_phys_cc_periodic(mmap, idir, idx + 1, gd) -
                             x_phys_cc_periodic(mmap, idir, idx - 1, gd));
}

#endif // PELELMEX_MESHMAP_EVALUATOR_H
