#ifndef AMREX_ML_ABECCECLAPLACIAN_H
#define AMREX_ML_ABECCECLAPLACIAN_H
#include <AMReX_Config.H>

#include <AMReX_MLCellABecCecLap.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad) - eta * (del dot c)) phi

class MLABecCecLaplacian : public MLCellABecCecLap
{
public:
  MLABecCecLaplacian() {}
  MLABecCecLaplacian(
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info = LPInfo(),
    const Vector<FabFactory<FArrayBox> const*>& a_factory = {},
    const int a_ncomp = 1);

  MLABecCecLaplacian(
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const Vector<iMultiFab const*>& a_overset_mask, // 1: unknown, 0: known
    const LPInfo& a_info = LPInfo(),
    const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

  virtual ~MLABecCecLaplacian();

  MLABecCecLaplacian(const MLABecCecLaplacian&) = delete;
  MLABecCecLaplacian(MLABecCecLaplacian&&) = delete;
  MLABecCecLaplacian& operator=(const MLABecCecLaplacian&) = delete;
  MLABecCecLaplacian& operator=(MLABecCecLaplacian&&) = delete;

  void define(
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info = LPInfo(),
    const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

  void define(
    const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const Vector<iMultiFab const*>& a_overset_mask,
    const LPInfo& a_info = LPInfo(),
    const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

  void setScalars(Real a, Real b, Real c) noexcept;
  void setRelaxation(Real omega) noexcept;
  void setACoeffs(int amrlev, const MultiFab& alpha);
  void setACoeffs(int amrlev, Real alpha);
  void
  setBCoeffs(int amrlev, const Array<MultiFab const*, AMREX_SPACEDIM>& beta);
  void setBCoeffs(int amrlev, Real beta);
  void setBCoeffs(int amrlev, Vector<Real> const& beta);
  void
  setCCoeffs(int amrlev, const Array<MultiFab const*, AMREX_SPACEDIM>& eta);

  void checkDiagonalDominance(int amrlev, int mglev);
  void getDiagonal(int amrlev, MultiFab& diag);
  void getNetFaceCoeff(MultiFab& coeff_d, int dir);

  virtual int getNComp() const override { return m_ncomp; }

  virtual bool needsUpdate() const override
  {
    return (m_needs_update || MLCellABecCecLap::needsUpdate());
  }
  virtual void update() override;

  virtual void prepareForSolve() override;
  virtual bool isSingular(int amrlev) const override
  {
    return m_is_singular[amrlev];
  }
  virtual bool isBottomSingular() const override { return m_is_singular[0]; }
  virtual void Fapply(int amrlev, int mglev, MultiFab& out, const MultiFab& in)
    const final override;
  virtual void Fsmooth(
    int amrlev,
    int mglev,
    MultiFab& sol,
    const MultiFab& rhs,
    int redblack) const final override;
  virtual void FFlux(
    int amrlev,
    const MFIter& mfi,
    const Array<FArrayBox*, AMREX_SPACEDIM>& flux,
    const FArrayBox& sol,
    Location /* loc */,
    const int face_only = 0) const final override;

  virtual void
  normalize(int amrlev, int mglev, MultiFab& mf) const final override;

  virtual Real getAScalar() const final override { return m_a_scalar; }
  virtual Real getBScalar() const final override { return m_b_scalar; }
  virtual Real getCScalar() const final override { return m_c_scalar; }
  virtual MultiFab const* getACoeffs(int amrlev, int mglev) const final override
  {
    return &(m_a_coeffs[amrlev][mglev]);
  }
  virtual Array<MultiFab const*, AMREX_SPACEDIM>
  getBCoeffs(int amrlev, int mglev) const final override
  {
    return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]);
  }
  virtual Array<MultiFab const*, AMREX_SPACEDIM>
  getCCoeffs(int amrlev, int mglev) const final override
  {
    return amrex::GetArrOfConstPtrs(m_c_coeffs[amrlev][mglev]);
  }

  virtual std::unique_ptr<MLLinOp>
  makeNLinOp(int /*grid_size*/) const final override
  {
    amrex::Abort("MLABecLaplacian::makeNLinOp: Not implemented");
    return std::unique_ptr<MLLinOp>{};
  }

  void averageDownCoeffsSameAmrLevel(
    int amrlev,
    Vector<MultiFab>& a,
    Vector<Array<MultiFab, AMREX_SPACEDIM>>& b,
    Vector<Array<MultiFab, AMREX_SPACEDIM>>& c);
  void averageDownCoeffs();
  void averageDownCoeffsToCoarseAmrLevel(int flev);

  void applyMetricTermsCoeffs();

  static void FFlux(
    Box const& box,
    Real const* dxinv,
    Real bscalar,
    Real cscalar,
    Array<FArrayBox const*, AMREX_SPACEDIM> const& bcoef,
    Array<FArrayBox const*, AMREX_SPACEDIM> const& ccoef,
    Array<FArrayBox*, AMREX_SPACEDIM> const& flux,
    FArrayBox const& sol,
    int face_only,
    int ncomp);

protected:
  bool m_needs_update = true;

  Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
  Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
  Real m_c_scalar = std::numeric_limits<Real>::quiet_NaN();
  Real m_omega = std::numeric_limits<Real>::quiet_NaN();
  Vector<Vector<MultiFab>> m_a_coeffs;
  Vector<Vector<Array<MultiFab, AMREX_SPACEDIM>>> m_b_coeffs;
  Vector<Vector<Array<MultiFab, AMREX_SPACEDIM>>> m_c_coeffs;

  Vector<int> m_is_singular;

  virtual bool supportRobinBC() const noexcept override { return false; }

private:
  int m_ncomp = 1;

  void define_abc_coeffs();
};

} // namespace amrex

#endif
