#ifndef PELELMEX_INTERIORSTRETCHMAP_H
#define PELELMEX_INTERIORSTRETCHMAP_H

#include <PeleLMeX_MeshMap.H>

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>

/**
 * \brief Interior-clustered mesh stretching with independent per-axis,
 *        per-side coarsening rates.
 *
 * Companion to \c TanhStretchMap: instead of clustering cells at the two
 * ends of an axis, this clusters them at an interior point and coarsens
 * outwards in both directions, at rates set independently by
 * \f$\beta_{\mathrm{lo}}\f$ (towards \c prob_lo) and
 * \f$\beta_{\mathrm{hi}}\f$ (towards \c prob_hi).
 *
 * \section design Construction
 *
 * The map is defined by its metric rather than the other way round, which
 * is what makes the two sides join smoothly.  With
 * \f[
 *   s = \frac{\xi - x_{\mathrm{lo}}}{L}\in[0,1],\qquad
 *   \eta = \begin{cases}
 *      (s_c - s)/s_c,      & s \le s_c \quad\text{(lower side)}\\[2pt]
 *      (s - s_c)/(1-s_c),  & s \ge s_c \quad\text{(upper side)}
 *   \end{cases}\;\in[0,1],
 * \f]
 * so that \f$\eta = 0\f$ at the clustering point and \f$\eta = 1\f$ at the
 * corresponding domain boundary, the stretch factor is \em prescribed as
 * \f[
 *   \mathrm{fac}_i(\xi) = \frac{dx_{\mathrm{phys}}}{d\xi}
 *       = \frac{\cosh^2\!\bigl(\beta_{\mathrm{side}}\,\eta\bigr)}
 *              {\mathcal{N}}.
 * \f]
 * Integrating, with \f$\;\mathcal{R}(\beta,\eta) \equiv
 * \int_0^\eta \cosh^2(\beta t)\,dt
 * = \dfrac{\eta}{2} + \dfrac{\sinh(2\beta\eta)}{4\beta}\;\f$ and
 * \f$\;\mathcal{G}(\beta)\equiv\mathcal{R}(\beta,1)\f$:
 * \f[
 *   x_{\mathrm{phys}}(\xi) = \begin{cases}
 *     x_c - L\,\dfrac{s_c}{\mathcal{N}}\;
 *           \mathcal{R}(\beta_{\mathrm{lo}},\eta), & s \le s_c,\\[10pt]
 *     x_c + L\,\dfrac{1-s_c}{\mathcal{N}}\;
 *           \mathcal{R}(\beta_{\mathrm{hi}},\eta), & s \ge s_c,
 *   \end{cases}
 *   \qquad
 *   \mathcal{N} = s_c\,\mathcal{G}(\beta_{\mathrm{lo}})
 *               + (1-s_c)\,\mathcal{G}(\beta_{\mathrm{hi}}).
 * \f]
 * The single \em global normalisation \f$\mathcal{N}\f$ is what preserves
 * the total extent: the physical extent equals the AMReX (Xi-space) extent
 * in every direction, exactly, for any pair of \f$\beta\f$'s, so as with
 * the other maps the mapping only redistributes cells within the domain
 * (mean \f$\mathrm{fac}_i = 1\f$).
 *
 * \section props Properties
 *
 * - The max/min cell-size ratio is \f$\cosh^2\beta_{\mathrm{lo}}\f$ below
 *   the clustering point and \f$\cosh^2\beta_{\mathrm{hi}}\f$ above it ---
 *   the same \f$\beta\f$ semantics as \c TanhStretchMap
 *   (\f$\beta = 2 \Rightarrow 14.2\,\times\f$,
 *   \f$\beta = 3 \Rightarrow 101\,\times\f$).
 * - Because \f$\cosh^2\f$ has a stationary minimum at \f$\eta = 0\f$, the
 *   metric is \f$C^1\f$ and the map \f$C^2\f$ across the join even when
 *   \f$\beta_{\mathrm{lo}} \neq \beta_{\mathrm{hi}}\f$: no jump in cell
 *   size and no kink in \c fac for the operators to trip over.  (Gluing
 *   two one-sided \c tanh maps at \f$\xi_c\f$ instead would leave a
 *   cell-size jump of
 *   \f$[\beta_{\mathrm{lo}}/\sinh 2\beta_{\mathrm{lo}}]/
 *      [\beta_{\mathrm{hi}}/\sinh 2\beta_{\mathrm{hi}}]\f$ there --- a
 *   factor of 3.4 for \f$(1.8, 2.6)\f$.)
 * - \f$\beta_{\mathrm{side}} = 0\f$ (or below 1e-8) gives uniform spacing
 *   on that side; both zero recovers the identity mapping along that axis
 *   exactly.
 * - The smallest cell on an axis carrying \c n cells is
 *   \f$\Delta x_{\min} = L/(n\,\mathcal{N})\f$.
 *
 * \section anchor Where the clustering point sits
 *
 * With independent \f$\beta\f$'s one cannot simultaneously have (i) the
 * total extent preserved, (ii) a continuous cell size at the join, and
 * (iii) \f$x_{\mathrm{phys}}(\xi_c) = \xi_c\f$.  This map keeps (i) and
 * (ii) and pins the \em physical location of the refinement, solving for
 * the Xi-space split in closed form.  With
 * \f$p = (x_c - x_{\mathrm{lo}})/L\f$:
 * \f[
 *   s_c = \frac{p\,\mathcal{G}(\beta_{\mathrm{hi}})}
 *              {(1-p)\,\mathcal{G}(\beta_{\mathrm{lo}})
 *               + p\,\mathcal{G}(\beta_{\mathrm{hi}})}.
 * \f]
 * \f$s_c\f$ is the fraction of the \em cells placed below \f$x_c\f$.  It
 * equals \f$p\f$ when \f$\beta_{\mathrm{lo}} = \beta_{\mathrm{hi}}\f$ and
 * otherwise shifts towards the \em less stretched side, which is the
 * physically correct behaviour: fast-growing cells need fewer of them to
 * span a given length.  Note that consequently \f$\xi_c \neq x_c\f$ in
 * general; \c splitXi() reports the split and \c centerXi() the
 * Xi-space coordinate of the clustering point.
 *
 * \warning The flip side of that shift is cell starvation.  If one
 * \f$\beta\f$ is much larger than the other, its side can be left with
 * only a handful of cells and the cell-to-cell growth becomes severe.  The
 * worst neighbouring-cell size ratio on a side carrying
 * \f$n_{\mathrm{side}}\f$ cells is
 * \f[
 *    g \simeq \exp\!\left(\frac{2\beta\tanh\beta}{n_{\mathrm{side}}}\right)
 *    \quad\Longrightarrow\quad
 *    n_{\mathrm{side}} \gtrsim \frac{2\beta\tanh\beta}{\ln g_{\max}} ,
 * \f]
 * i.e. \f$\beta = 3\f$ at a 10% growth limit needs \f$\gtrsim 63\f$ cells
 * on that side.  The constructor reads \c amr.n_cell and prints a warning
 * when the estimate exceeds \c InteriorStretchMap.max_growth (default
 * 1.15) or when a side is left with fewer than 4 cells on the base level.
 *
 * ParmParse block:
 * \code
 *   # coarsening rates, one value per AMREX_SPACEDIM axis, both >= 0
 *   InteriorStretchMap.beta_lo     = 1.8 0.0
 *   InteriorStretchMap.beta_hi     = 2.6 0.0
 *
 *   # clustering point: EITHER physical coordinates ...
 *   InteriorStretchMap.center      = 0.03 0.0
 *   # ... OR fractions of the axis in (0,1); default 0.5 (mid-axis)
 *   InteriorStretchMap.center_frac = 0.3 0.5
 *
 *   # growth limit used only for the start-up warning; default 1.15
 *   InteriorStretchMap.max_growth  = 1.15
 * \endcode
 * \c center wins if both are given.  On an axis where both betas vanish
 * the clustering point is irrelevant (the map is the identity there).
 */
class InteriorStretchMap : public MeshMap
{
public:
  InteriorStretchMap();
  ~InteriorStretchMap() override = default;
  InteriorStretchMap(const InteriorStretchMap&) = delete;
  InteriorStretchMap& operator=(const InteriorStretchMap&) = delete;
  InteriorStretchMap(InteriorStretchMap&&) = delete;
  InteriorStretchMap& operator=(InteriorStretchMap&&) = delete;

  void create_map(int lev, const amrex::Geometry& geom) override;

  [[nodiscard]] MeshMapEvaluator make_evaluator() const override;

  void fill_nodal_displacement(
    int lev,
    const amrex::Geometry& geom,
    amrex::MultiFab& disp_nd) const override;

  // Accessors (useful for drivers / diagnostic code that wants to replay
  // the mapping formula, e.g. when writing mapping-aware initial
  // conditions).

  /// Coarsening rate towards \c prob_lo on axis \p dir.
  [[nodiscard]] amrex::Real beta_lo(int dir) const noexcept
  {
    return m_beta_lo[dir];
  }
  /// Coarsening rate towards \c prob_hi on axis \p dir.
  [[nodiscard]] amrex::Real beta_hi(int dir) const noexcept
  {
    return m_beta_hi[dir];
  }
  /// Clustering point on axis \p dir as a fraction of the axis, in (0,1).
  [[nodiscard]] amrex::Real centerFrac(int dir) const noexcept
  {
    return m_center_frac[dir];
  }
  /// Xi-space split on axis \p dir: the fraction of cells below the
  /// clustering point.
  [[nodiscard]] amrex::Real splitXi(int dir) const noexcept
  {
    return m_s_c[dir];
  }
  /// Normalisation \f$\mathcal{N}\f$ on axis \p dir; \c fac_min = 1/N.
  [[nodiscard]] amrex::Real normalization(int dir) const noexcept
  {
    return m_N[dir];
  }
  /// Physical coordinate of the clustering point on axis \p dir.
  [[nodiscard]] amrex::Real center(int dir) const noexcept
  {
    return m_prob_lo[dir] + m_Lxi[dir] * m_center_frac[dir];
  }
  /// Xi-space coordinate of the clustering point on axis \p dir; equals
  /// \c center(dir) only when \c beta_lo == \c beta_hi.
  [[nodiscard]] amrex::Real centerXi(int dir) const noexcept
  {
    return m_prob_lo[dir] + m_Lxi[dir] * m_s_c[dir];
  }

private:
  amrex::Real m_beta_lo[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};
  amrex::Real m_beta_hi[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};

  // Clustering point as a fraction of the axis, in (0,1).
  amrex::Real m_center_frac[AMREX_SPACEDIM] = {AMREX_D_DECL(0.5, 0.5, 0.5)};

  // Derived at construction from (beta_lo, beta_hi, center_frac).
  amrex::Real m_s_c[AMREX_SPACEDIM] = {AMREX_D_DECL(0.5, 0.5, 0.5)};
  amrex::Real m_N[AMREX_SPACEDIM] = {AMREX_D_DECL(1.0, 1.0, 1.0)};

  // Domain, cached from the geometry ParmParse block so that the
  // accessors and the start-up report can be expressed in physical units
  // without a Geometry in hand.
  amrex::Real m_prob_lo[AMREX_SPACEDIM] = {AMREX_D_DECL(0.0, 0.0, 0.0)};
  amrex::Real m_Lxi[AMREX_SPACEDIM] = {AMREX_D_DECL(1.0, 1.0, 1.0)};
};

#endif /* PELELMEX_INTERIORSTRETCHMAP_H */
