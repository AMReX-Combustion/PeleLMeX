#ifndef PELELMEX_MESH_MAPPED_CELL_CONS_INTERP_H
#define PELELMEX_MESH_MAPPED_CELL_CONS_INTERP_H

#include <AMReX_MFInterpolater.H>
#include <AMReX_MultiFab.H>

/**
 * \brief Volume-conservative cell-centered velocity interpolator for
 *        the mesh-mapping path.
 *
 * Wraps \c amrex::mf_cell_cons_interp by pre-weighting the coarse
 * velocity into Xi-space (\f$u^{xi}_i = v^{phys}_i \cdot \mathrm{detJ}
 * /\mathrm{fac}_i\f$), running the standard cell-conservative linear
 * interp on that representation (which is volume-conservative since
 * Xi-space cells are uniform), and un-weighting on the fine side.
 *
 * Assumes the data is per-direction velocity with exactly
 * \c AMREX_SPACEDIM components.  Used by \c PeleLM::getInterpolator()
 * in place of \c amrex::mf_cell_cons_interp when mesh mapping is
 * active.  Lifetime managed by \c PeleLM via the \c m_mapped_interps
 * vector (rebuilt on every level lifecycle event) -- AMReX's
 * \c FabArrayBase::TheFPinfo cache holds a raw pointer to the
 * interpolator and requires it to outlive a single fillpatch call.
 */
class MeshMappedCellConsInterp final : public amrex::MFInterpolater
{
public:
  MeshMappedCellConsInterp(
    amrex::MultiFab const* a_detJ_crse,
    amrex::MultiFab const* a_fac_crse,
    amrex::MultiFab const* a_detJ_fine,
    amrex::MultiFab const* a_fac_fine)
    : m_detJ_crse(a_detJ_crse),
      m_fac_crse(a_fac_crse),
      m_detJ_fine(a_detJ_fine),
      m_fac_fine(a_fac_fine)
  {
  }

  amrex::Box CoarseBox(amrex::Box const& fine, int ratio) override
  {
    return amrex::mf_cell_cons_interp.CoarseBox(fine, ratio);
  }
  amrex::Box
  CoarseBox(amrex::Box const& fine, amrex::IntVect const& ratio) override
  {
    return amrex::mf_cell_cons_interp.CoarseBox(fine, ratio);
  }

  void interp(
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
    int bcscomp) override;

private:
  // Non-owning; lifetime managed by PeleLM (rebuilt alongside MeshMap
  // metric MultiFabs on every level lifecycle event).
  amrex::MultiFab const* m_detJ_crse{nullptr};
  amrex::MultiFab const* m_fac_crse{nullptr};
  amrex::MultiFab const* m_detJ_fine{nullptr};
  amrex::MultiFab const* m_fac_fine{nullptr};
};

#endif // PELELMEX_MESH_MAPPED_CELL_CONS_INTERP_H
