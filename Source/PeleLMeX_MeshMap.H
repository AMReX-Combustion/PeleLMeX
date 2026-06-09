#ifndef PELELMEX_MESHMAP_H
#define PELELMEX_MESHMAP_H

#include <memory>
#include <string>

#include <AMReX_Array.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_FabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_SPACE.H>
#include <AMReX_Vector.H>

#include <PeleLMeX_MeshMapEvaluator.H>

/**
 * \brief Abstract mesh-mapping base class.
 *
 * The MeshMap object owns per-level storage for the metric fields associated
 * with a single-block, diagonal, axis-aligned mesh mapping
 *     X = Phi(Xi) = (x(chi), y(eta), z(zeta))
 * on top of the uniform AMReX mesh.  Faces remain axis-aligned; the only
 * new quantities are the per-direction stretch factors fac_i and the
 * Jacobian detJ = prod_i fac_i.
 *
 * Each level stores, at each of CC / ND / XF / YF / ZF locations:
 *   - fac  : AMREX_SPACEDIM-component field giving per-direction stretch
 *   - detJ : 1-component field giving the Jacobian determinant
 *
 * The concrete class implements create_map(lev, geom) which fills these
 * metric fields from the geometry.
 *
 * Mirrors the design of amr-wind/core/MeshMap.H (PR #545 in Exawind/amr-wind)
 * but stores MultiFabs directly on the object rather than going through a
 * FieldRepo, since PeleLMeX uses raw AMReX MultiFabs for its state.
 */
class MeshMap
{
public:
  MeshMap() = default;
  virtual ~MeshMap() = default;
  MeshMap(const MeshMap&) = delete;
  MeshMap& operator=(const MeshMap&) = delete;
  MeshMap(MeshMap&&) = delete;
  MeshMap& operator=(MeshMap&&) = delete;

  /**
   * \brief Resize the per-level storage vectors to hold \p nlev entries.
   *
   * Called by PeleLM when the level count changes (after MakeNewLevel,
   * ClearLevel, or regrid).  Individual levels are populated via define().
   */
  void resize(int nlev);

  /**
   * \brief Clear the MultiFabs stored at level \p lev.
   *
   * Equivalent to amrex::MultiFab::clear() on every per-location field.
   * Invoked from PeleLM::ClearLevel().
   */
  void clear_level(int lev);

  /**
   * \brief Define (allocate) all metric MultiFabs on level \p lev.
   *
   * \p nghost is the number of ghost cells to request; matches the
   * maximum ghost-cell requirement of the consumers (projectors,
   * advection, diffusion).
   */
  void define(
    int lev,
    const amrex::BoxArray& ba,
    const amrex::DistributionMapping& dm,
    const amrex::FabFactory<amrex::FArrayBox>& factory,
    int nghost);

  /**
   * \brief Fill the metric fields on level \p lev from \p geom.
   *
   * Must be implemented by every concrete derived class.
   */
  virtual void create_map(int lev, const amrex::Geometry& geom) = 0;

  /**
   * \brief Fill a nodal MultiFab with the per-direction displacement
   *        from the Xi-space node position to the physical-space node
   *        position.
   *
   * Used when writing plotfiles with mesh-mapping metadata so that
   * ParaView / VisIt's AMReX reader can render the solution on the
   * curvilinear physical grid (see PeleLMeX_Plot.cpp's
   * writePlotfileWithMapping).  The displacement convention is
   *     disp_i(node) = x_phys_i(node) - x_xi_i(node)
   * The default implementation assumes a piecewise-constant stretch
   * factor (diagonal, axis-aligned, single-block) and reads the
   * per-direction factors out of the already-filled nodal fac MultiFab:
   *     disp_i(node) = (fac_i(node) - 1) * (x_xi_i(node) - prob_lo_i)
   * which is exact for ConstantMap.  Derived classes with spatially
   * varying fac must override this method with the correct integral
   *     disp_i(node) = int_{prob_lo_i}^{x_xi_i(node)} (fac_i(xi) - 1) dxi_i
   *
   * \p disp_nd must be a nodal MultiFab with at least AMREX_SPACEDIM
   * components.  Only the first AMREX_SPACEDIM components are written.
   */
  virtual void fill_nodal_displacement(
    int lev, const amrex::Geometry& geom, amrex::MultiFab& disp_nd) const;

  /**
   * \brief Return a device-callable evaluator capable of producing the
   *        physical-space coordinate of any cell centre or face from
   *        its integer index plus the corresponding GeometryData.
   *
   * The base implementation returns the identity (\c Kind::Identity);
   * each concrete map populates a richer payload.  Held by value by
   * the PeleLM driver and captured by value into device lambdas, this
   * is the recommended way for user \c initdata_mapped /
   * \c bcnormal code to obtain physical positions in a mapping-aware
   * but mapping-agnostic fashion.
   */
  [[nodiscard]] virtual MeshMapEvaluator make_evaluator() const
  {
    return MeshMapEvaluator{};
  }

  /** \return the number of active levels currently allocated. */
  [[nodiscard]] int num_levels() const
  {
    return static_cast<int>(m_fac_cc.size());
  }

  // ---- Cell-centered accessors -------------------------------------------
  amrex::MultiFab& fac_cc(int lev) { return m_fac_cc[lev]; }
  amrex::MultiFab& detJ_cc(int lev) { return m_detJ_cc[lev]; }
  [[nodiscard]] const amrex::MultiFab& fac_cc(int lev) const
  {
    return m_fac_cc[lev];
  }
  [[nodiscard]] const amrex::MultiFab& detJ_cc(int lev) const
  {
    return m_detJ_cc[lev];
  }

  // ---- Nodal accessors ---------------------------------------------------
  amrex::MultiFab& fac_nd(int lev) { return m_fac_nd[lev]; }
  amrex::MultiFab& detJ_nd(int lev) { return m_detJ_nd[lev]; }
  [[nodiscard]] const amrex::MultiFab& fac_nd(int lev) const
  {
    return m_fac_nd[lev];
  }
  [[nodiscard]] const amrex::MultiFab& detJ_nd(int lev) const
  {
    return m_detJ_nd[lev];
  }

  // ---- Face-centered accessors (per direction) --------------------------
  amrex::MultiFab& fac_fc(int lev, int idim) { return m_fac_fc[lev][idim]; }
  amrex::MultiFab& detJ_fc(int lev, int idim) { return m_detJ_fc[lev][idim]; }
  [[nodiscard]] const amrex::MultiFab& fac_fc(int lev, int idim) const
  {
    return m_fac_fc[lev][idim];
  }
  [[nodiscard]] const amrex::MultiFab& detJ_fc(int lev, int idim) const
  {
    return m_detJ_fc[lev][idim];
  }

  /// Array accessor: returns references to the per-direction face MFs on a
  /// level.
  amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> fac_fc(int lev)
  {
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> r;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      r[d] = &m_fac_fc[lev][d];
    }
    return r;
  }
  amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> detJ_fc(int lev)
  {
    amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> r;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
      r[d] = &m_detJ_fc[lev][d];
    }
    return r;
  }

  /**
   * \brief Factory method: build a concrete MeshMap by name.
   *
   * Recognised names:
   *   - "ConstantMap"    : uniform per-direction stretching
   *   - "ExpStretchMap"  : one-sided (single-axis) exponential mesh stretching
   *   - "TanhStretchMap" : two-sided tanh-based mesh stretching
   *
   * Aborts with an informative message if the name is not recognised.
   */
  static std::unique_ptr<MeshMap> create(const std::string& name);

protected:
  amrex::Vector<amrex::MultiFab> m_fac_cc;
  amrex::Vector<amrex::MultiFab> m_detJ_cc;

  amrex::Vector<amrex::MultiFab> m_fac_nd;
  amrex::Vector<amrex::MultiFab> m_detJ_nd;

  amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> m_fac_fc;
  amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> m_detJ_fc;
};

#endif /* PELELMEX_MESHMAP_H */
