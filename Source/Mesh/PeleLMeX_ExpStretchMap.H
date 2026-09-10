#ifndef PELELMEX_EXPSTRETCHMAP_H
#define PELELMEX_EXPSTRETCHMAP_H

#include <PeleLMeX_MeshMap.H>

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

/**
 * \brief One-sided exponential mesh stretching along a single axis.
 *
 * The stretched direction idir is mapped by
 * \f[
 *    x_{\mathrm{phys}}(\eta) = x_{\mathrm{lo}} + L\,
 *      \frac{e^{\beta\eta} - 1}{e^\beta - 1} \quad\text{(wall at lo)}
 * \f]
 * or, mirrored, for wall at hi, with the normalized coordinate
 * \f$\eta \in [0,1]\f$ measuring distance from the wall end.  The AMReX
 * (Xi-space) domain spans the same physical extent as \p prob_lo /
 * \p prob_hi, so the total length is preserved and the mapping just
 * redistributes cells within it (mean \f$\mathrm{fac}_i = 1\f$).  The
 * stretch factor is
 * \f[
 *    \mathrm{fac}(\eta) = \frac{\beta\,e^{\beta\eta}}{e^\beta - 1}
 * \f]
 * with \f$\mathrm{fac}(0) = \beta/(e^\beta - 1)\f$ at the wall (smallest
 * cells) and \f$\mathrm{fac}(1) = \beta e^\beta /(e^\beta - 1)\f$ far
 * from the wall; the ratio of largest to smallest physical cell is
 * exactly \f$e^\beta\f$.  The non-stretched directions carry
 * \f$\mathrm{fac} = 1\f$ everywhere.
 *
 * ParmParse block:
 * \code
 *   ExpStretchMap.direction = 1          # axis (0,1,2); default 1 (=y)
 *   ExpStretchMap.wall_end  = lo         # "lo" or "hi"; default "lo"
 *   ExpStretchMap.beta      = 2.0        # stretching strength >= 0
 * \endcode
 *
 * \f$\beta = 0\f$ (or < 1e-8) recovers the identity mapping.
 */
class ExpStretchMap : public MeshMap
{
public:
  ExpStretchMap();
  ~ExpStretchMap() override = default;
  ExpStretchMap(const ExpStretchMap&) = delete;
  ExpStretchMap& operator=(const ExpStretchMap&) = delete;
  ExpStretchMap(ExpStretchMap&&) = delete;
  ExpStretchMap& operator=(ExpStretchMap&&) = delete;

  void create_map(int lev, const amrex::Geometry& geom) override;

  [[nodiscard]] MeshMapEvaluator make_evaluator() const override;

  void fill_nodal_displacement(
    int lev,
    const amrex::Geometry& geom,
    amrex::MultiFab& disp_nd) const override;

  // Accessors (useful for drivers / diagnostic code that wants to replay
  // the mapping formula, e.g. when writing mapping-aware initial
  // conditions).
  [[nodiscard]] int direction() const noexcept { return m_dir; }
  [[nodiscard]] int wall_end() const noexcept { return m_wall_end; }
  [[nodiscard]] amrex::Real beta() const noexcept { return m_beta; }

private:
  int m_dir{1};      // stretched direction (0,1,2); default y
  int m_wall_end{0}; // 0 = lo, 1 = hi
  amrex::Real m_beta{1.0};
};

#endif /* PELELMEX_EXPSTRETCHMAP_H */
