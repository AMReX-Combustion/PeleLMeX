#ifndef PELELMEX_TANHSTRETCHMAP_H
#define PELELMEX_TANHSTRETCHMAP_H

#include <PeleLMeX_MeshMap.H>

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

/**
 * \brief Two-sided tanh-based mesh stretching with independent per-axis
 *        clustering.
 *
 * Each axis i is stretched by the length-preserving mapping
 * \f[
 *    x_{\mathrm{phys}}(\xi) = x_{\mathrm{lo}} + \frac{L}{2}\left(
 *       1 + \frac{\tanh\bigl(\beta_i(2\xi_n - 1)\bigr)}{\tanh\beta_i}
 *    \right),\qquad \xi_n = \frac{\xi - x_{\mathrm{lo}}}{L}\in[0,1].
 * \f]
 * The corresponding stretch factor is
 * \f[
 *    \mathrm{fac}_i(\xi) = \frac{dx_{\mathrm{phys}}}{d\xi}
 *        = \frac{\beta_i\,\mathrm{sech}^2\!\bigl(\beta_i(2\xi_n - 1)\bigr)}
 *               {\tanh\beta_i}.
 * \f]
 * For \f$\beta_i > 0\f$ this clusters cells at \em both ends of the
 * axis with a max/min cell-size ratio of \f$\cosh^2\beta_i\f$
 * (e.g. \f$\beta = 2 \Rightarrow 14.2\,\times\f$,
 * \f$\beta = 3 \Rightarrow 101\,\times\f$).  Setting \f$\beta_i = 0\f$
 * (or below 1e-8) recovers the identity mapping along that axis.
 *
 * Total physical extent equals the AMReX (Xi-space) extent in every
 * direction; the mapping only redistributes cells within the domain.
 *
 * ParmParse block:
 * \code
 *   TanhStretchMap.beta = 2.0 2.0 0.0   # one float per AMREX_SPACEDIM
 * \endcode
 * Each entry is non-negative; entries below 1e-8 are treated as
 * uniform.
 */
class TanhStretchMap : public MeshMap
{
public:
  TanhStretchMap();
  ~TanhStretchMap() override = default;
  TanhStretchMap(const TanhStretchMap&) = delete;
  TanhStretchMap& operator=(const TanhStretchMap&) = delete;
  TanhStretchMap(TanhStretchMap&&) = delete;
  TanhStretchMap& operator=(TanhStretchMap&&) = delete;

  void create_map(int lev, const amrex::Geometry& geom) override;

  [[nodiscard]] MeshMapEvaluator make_evaluator() const override;

  void fill_nodal_displacement(
    int lev,
    const amrex::Geometry& geom,
    amrex::MultiFab& disp_nd) const override;

  // Accessors (useful for drivers / diagnostic code that wants to
  // replay the mapping formula, e.g. when writing mapping-aware
  // initial conditions).
  [[nodiscard]] amrex::Real beta(int dir) const noexcept { return m_beta[dir]; }

private:
  amrex::Real m_beta[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};
};

#endif /* PELELMEX_TANHSTRETCHMAP_H */
