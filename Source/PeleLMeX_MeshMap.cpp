#include <PeleLMeX_MeshMap.H>
#include <PeleLMeX_ConstantMap.H>

#include <AMReX.H>
#include <AMReX_IntVect.H>
#include <AMReX_MultiFab.H>

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
  m_fac_cc[lev].define(
    ba, dm, AMREX_SPACEDIM, ng, amrex::MFInfo(), factory);
  m_detJ_cc[lev].define(ba, dm, 1, ng, amrex::MFInfo(), factory);

  // Nodal.
  const amrex::BoxArray nba = amrex::convert(ba, amrex::IntVect::TheNodeVector());
  m_fac_nd[lev].define(
    nba, dm, AMREX_SPACEDIM, ng, amrex::MFInfo(), factory);
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
  amrex::Abort(
    "MeshMap::create(): unrecognised mesh-mapping name '" + name +
    "'.  Supported: ConstantMap.");
  return nullptr;
}
