#ifndef AMREX_ML_CELL_ABECCECLAP_H
#define AMREX_ML_CELL_ABECCECLAP_H
#include <AMReX_Config.H>

#include <AMReX_MLCellLinOp.H>

namespace amrex {

class MLCellABecCecLap : public MLCellLinOpT<MultiFab>
{
public:
  MLCellABecCecLap();
  virtual ~MLCellABecCecLap();

  MLCellABecCecLap(const MLCellABecCecLap&) = delete;
  MLCellABecCecLap(MLCellABecCecLap&&) = delete;
  MLCellABecCecLap& operator=(const MLCellABecCecLap&) = delete;
  MLCellABecCecLap& operator=(MLCellABecCecLap&&) = delete;

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

  iMultiFab const* getOversetMask(int amrlev, int mglev) const
  {
    return m_overset_mask[amrlev][mglev].get();
  }

  virtual bool needsUpdate() const override
  {
    return MLCellLinOp::needsUpdate();
  }
  virtual void update() override;

  virtual void prepareForSolve() override;

  virtual void getFluxes(
    const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& a_flux,
    const Vector<MultiFab*>& a_sol,
    Location a_loc) const final override;
  virtual void getFluxes(
    const Vector<MultiFab*>& /*a_flux*/,
    const Vector<MultiFab*>& /*a_sol*/) const final override
  {
    amrex::Abort("MLCellABecLap::getFluxes: How did we get here?");
  }

  virtual Real getAScalar() const = 0;
  virtual Real getBScalar() const = 0;
  virtual Real getCScalar() const = 0;
  virtual MultiFab const* getACoeffs(int amrlev, int mglev) const = 0;
  virtual Array<MultiFab const*, AMREX_SPACEDIM>
  getBCoeffs(int amrlev, int mglev) const = 0;
  virtual Array<MultiFab const*, AMREX_SPACEDIM>
  getCCoeffs(int amrlev, int mglev) const = 0;

  virtual void
  applyInhomogNeumannTerm(int amrlev, MultiFab& rhs) const final override;

  virtual void applyOverset(int amlev, MultiFab& rhs) const override;
#ifdef AMREX_USE_HYPRE
  virtual std::unique_ptr<Hypre>
  makeHypre(Hypre::Interface hypre_interface) const override;
#endif

#ifdef AMREX_USE_PETSC
  virtual std::unique_ptr<PETScABecLap> makePETSc() const override;
#endif

protected:
  Vector<Vector<std::unique_ptr<iMultiFab>>> m_overset_mask;

  LPInfo m_lpinfo_arg;

  virtual bool supportInhomogNeumannBC() const noexcept override
  {
    return false;
  }
};

} // namespace amrex

#endif
