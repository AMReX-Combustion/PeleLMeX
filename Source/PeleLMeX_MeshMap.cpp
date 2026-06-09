#include <PeleLMeX_MeshMap.H>
#include <PeleLMeX_ConstantMap.H>
#include <PeleLMeX_ExpStretchMap.H>
#include <PeleLMeX_TanhStretchMap.H>

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_Gpu.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>

#include <memory>
#include <string>

void
MeshMap::resize(int nlev)
{
  m_fac_cc.resize(nlev);
  m_detJ_cc.resize(nlev);
  m_fac_nd.resize(nlev);
  m_detJ_nd.resize(nlev);
  m_fac_fc.resize(nlev);
  m_detJ_fc.resize(nlev);
}

void
MeshMap::clear_level(int lev)
{
  m_fac_cc[lev].clear();
  m_detJ_cc[lev].clear();
  m_fac_nd[lev].clear();
  m_detJ_nd[lev].clear();
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    m_fac_fc[lev][d].clear();
    m_detJ_fc[lev][d].clear();
  }
}

void
MeshMap::define(
  int lev,
  const amrex::BoxArray& ba,
  const amrex::DistributionMapping& dm,
  const amrex::FabFactory<amrex::FArrayBox>& factory,
  int nghost)
{
  AMREX_ASSERT(lev >= 0 && lev < num_levels());
  const amrex::IntVect ng(nghost);

  // Cell-centered: fac has AMREX_SPACEDIM components, detJ has 1.
  m_fac_cc[lev].define(ba, dm, AMREX_SPACEDIM, ng, amrex::MFInfo(), factory);
  m_detJ_cc[lev].define(ba, dm, 1, ng, amrex::MFInfo(), factory);

  // Nodal.
  const amrex::BoxArray nba =
    amrex::convert(ba, amrex::IntVect::TheNodeVector());
  m_fac_nd[lev].define(nba, dm, AMREX_SPACEDIM, ng, amrex::MFInfo(), factory);
  m_detJ_nd[lev].define(nba, dm, 1, ng, amrex::MFInfo(), factory);

  // Face-centered, per direction.  Matches the amr-wind convention of
  // storing all AMREX_SPACEDIM stretch components on every face even
  // though only one is strictly needed per face direction --- downstream
  // code consumes fac_fc[lev][idim](i,j,k,n) with n indexing direction.
  for (int d = 0; d < AMREX_SPACEDIM; ++d) {
    const amrex::BoxArray fba =
      amrex::convert(ba, amrex::IntVect::TheDimensionVector(d));
    m_fac_fc[lev][d].define(
      fba, dm, AMREX_SPACEDIM, ng, amrex::MFInfo(), factory);
    m_detJ_fc[lev][d].define(fba, dm, 1, ng, amrex::MFInfo(), factory);
  }
}

std::unique_ptr<MeshMap>
MeshMap::create(const std::string& name)
{
  if (name == "ConstantMap") {
    return std::make_unique<ConstantMap>();
  }
  if (name == "ExpStretchMap") {
    return std::make_unique<ExpStretchMap>();
  }
  if (name == "TanhStretchMap") {
    return std::make_unique<TanhStretchMap>();
  }
  amrex::Abort(
    "MeshMap::create(): unrecognised mesh-mapping name '" + name +
    "'.  Supported: ConstantMap, ExpStretchMap, TanhStretchMap.");
  return nullptr;
}

void
MeshMap::fill_nodal_displacement(
  int lev, const amrex::Geometry& geom, amrex::MultiFab& disp_nd) const
{
  AMREX_ASSERT(lev >= 0 && lev < num_levels());
  AMREX_ASSERT(disp_nd.nComp() >= AMREX_SPACEDIM);
  AMREX_ASSERT_WITH_MESSAGE(
    disp_nd.boxArray().ixType().nodeCentered(),
    "MeshMap::fill_nodal_displacement: output MultiFab must be nodal");

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> plo{
    AMREX_D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2))};
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx_xi{
    AMREX_D_DECL(geom.CellSize(0), geom.CellSize(1), geom.CellSize(2))};

  // disp_i(node) = (fac_i(node) - 1) * (x_xi_i(node) - prob_lo_i)
  //              = (fac_i(node) - 1) * i * dx_xi_i     (for node index i)
  // which is exact for ConstantMap (piecewise-constant fac).
  const auto& fac_ma = m_fac_nd[lev].const_arrays();
  const auto& disp_ma = disp_nd.arrays();
  amrex::ParallelFor(
    disp_nd, amrex::IntVect(0),
    [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
      auto f = fac_ma[box_no];
      auto d = disp_ma[box_no];
      AMREX_D_TERM(d(i, j, k, 0) = (f(i, j, k, 0) - amrex::Real(1.0)) *
                                   static_cast<amrex::Real>(i) * dx_xi[0];
                   , d(i, j, k, 1) = (f(i, j, k, 1) - amrex::Real(1.0)) *
                                     static_cast<amrex::Real>(j) * dx_xi[1];
                   , d(i, j, k, 2) = (f(i, j, k, 2) - amrex::Real(1.0)) *
                                     static_cast<amrex::Real>(k) * dx_xi[2];);
      amrex::ignore_unused(plo); // not needed in this formulation
    });
  amrex::Gpu::streamSynchronize();
}
