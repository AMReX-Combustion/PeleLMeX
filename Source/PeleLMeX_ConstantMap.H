#ifndef PELELMEX_CONSTANTMAP_H
#define PELELMEX_CONSTANTMAP_H

#include <PeleLMeX_MeshMap.H>

#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

/**
 * \brief Spatially uniform (per-direction) stretching mesh map.
 *
 * Fills the metric fields with a constant per-direction factor
 * \f$\mathrm{fac}_i\f$ specified via ParmParse:
 *
 * \code
 *   ConstantMap.scaling_factor = fac_x fac_y fac_z
 * \endcode
 *
 * Jacobian: \f$\det J = \prod_i \mathrm{fac}_i\f$.
 *
 * This is primarily useful as a sanity check: a uniform scaling should
 * reproduce the unmapped solution exactly in the incompressible limit
 * (with correctly scaled domain), providing an easy regression test for
 * the mesh-mapping machinery.
 */
class ConstantMap : public MeshMap
{
public:
  ConstantMap();
  ~ConstantMap() override = default;
  ConstantMap(const ConstantMap&) = delete;
  ConstantMap& operator=(const ConstantMap&) = delete;
  ConstantMap(ConstantMap&&) = delete;
  ConstantMap& operator=(ConstantMap&&) = delete;

  void create_map(int lev, const amrex::Geometry& geom) override;

  [[nodiscard]] MeshMapEvaluator make_evaluator() const override;

  void fill_cc(int lev);
  void fill_nd(int lev);
  void fill_fc(int lev);

private:
  amrex::Vector<amrex::Real> m_fac{AMREX_D_DECL(1.0, 1.0, 1.0)};
};

#endif /* PELELMEX_CONSTANTMAP_H */
