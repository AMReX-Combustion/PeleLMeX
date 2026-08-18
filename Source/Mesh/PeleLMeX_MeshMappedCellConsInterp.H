#ifndef PELELMEX_MESH_MAPPED_CELL_CONS_INTERP_H
#define PELELMEX_MESH_MAPPED_CELL_CONS_INTERP_H

#include <AMReX_MFInterpolater.H>
#include <AMReX_MultiFab.H>

/**
 * \brief Volume-conservative cell-centered state interpolator for the
 *        mesh-mapping path.
 *
 * Wraps \c amrex::mf_cell_cons_interp: stage the coarse data in Xi-space
 * (v*w), run the standard cell-cons interp (volume-conservative on the
 * uniform Xi mesh), then un-weight the fine side (v_xi/w).  Per-component
 * weight w: velocity [vel_start, vel_start+vel_ncomp) -> detJ/fac_i;
 * conserved densities [dens_start, dens_end] (rho, rhoY, rhoH) -> detJ;
 * everything else (TEMP, RHORT, ...) -> 1.
 *
 * Assumes the fillpatch starts at state component 0 (ccomp == fcomp == 0),
 * as PeleLM::fillpatch_state / fillcoarsepatch_state invoke it.  Components
 * left at identity weight (TEMP, RHORT) are reconciled by the EOS recompute
 * after prolongation in the regrid path.
 *
 * Used by PeleLM::getInterpolator() when mesh mapping is active.  Lifetime
 * is managed via PeleLM's m_mapped_interps (rebuilt on level-lifecycle
 * events) because AMReX's FabArrayBase::TheFPinfo cache holds a raw pointer
 * that must outlive a fillpatch call.
 */
class MeshMappedCellConsInterp final : public amrex::MFInterpolater
{
public:
  MeshMappedCellConsInterp(
    amrex::MultiFab const* a_detJ_crse,
    amrex::MultiFab const* a_fac_crse,
    amrex::MultiFab const* a_detJ_fine,
    amrex::MultiFab const* a_fac_fine,
    int a_vel_start,
    int a_vel_ncomp,
    int a_dens_start,
    int a_dens_end)
    : m_detJ_crse(a_detJ_crse),
      m_fac_crse(a_fac_crse),
      m_detJ_fine(a_detJ_fine),
      m_fac_fine(a_fac_fine),
      m_vel_start(a_vel_start),
      m_vel_ncomp(a_vel_ncomp),
      m_dens_start(a_dens_start),
      m_dens_end(a_dens_end)
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

  // Component layout for the per-component weight (see class doc); set
  // by PeleLM from VELX, AMREX_SPACEDIM, DENSITY, RHOH.
  int m_vel_start{0};
  int m_vel_ncomp{0};
  int m_dens_start{-1};
  int m_dens_end{-2}; // empty range by default (dens_end < dens_start)
};

#endif // PELELMEX_MESH_MAPPED_CELL_CONS_INTERP_H
