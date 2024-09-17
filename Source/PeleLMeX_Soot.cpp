
#ifdef PELE_USE_SOOT
#include <PeleLMeX.H>
#include <PeleLMeX_Derive.H>
#include "SootModel.H"

using namespace amrex;

void
PeleLM::setSootIndx()
{
  const int ndim = AMREX_SPACEDIM;
  SootComps sootComps;
  // U/Q state variables indices
  sootComps.qRhoIndx = DENSITY - ndim;
  sootComps.qSpecIndx = FIRSTSPEC - ndim;
  sootComps.qTempIndx = TEMP - ndim;
  sootComps.qSootIndx = FIRSTSOOT - ndim;
  // Source variables indices
  sootComps.rhoIndx = DENSITY - ndim;
  sootComps.specIndx = FIRSTSPEC - ndim;
  sootComps.engIndx = RHOH - ndim;
  sootComps.sootIndx = FIRSTSOOT - ndim;
  soot_model->setIndices(sootComps);
}

void
PeleLM::cleanupSootModel()
{
  soot_model->cleanup();
  delete soot_model;
}

void
PeleLM::computeSootSource(const PeleLM::TimeStamp& a_timestamp, const Real a_dt)
{
  bool pres_term = false; // Do not include change in pressure in energy
  for (int lev = 0; lev <= finest_level; lev++) {
    auto* ldata_p = getLevelDataPtr(lev, a_timestamp);
    Real time = getTime(lev, a_timestamp);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*(m_extSource[lev]), TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      Box const& gbx = mfi.growntilebox();
      auto const& state_arr = ldata_p->state.const_array(mfi, DENSITY);
      auto const& mu = ldata_p->visc_cc.const_array(mfi, 0);
      auto const& source_arr = m_extSource[lev]->array(mfi, DENSITY);
      soot_model->computeSootSourceTerm(
        gbx, state_arr, mu, source_arr, time, a_dt, pres_term);
    }
  }
}

void
PeleLM::clipSootMoments()
{
  for (int lev = 0; lev <= finest_level; lev++) {
    auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(ldata_p->state, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      Box const& gbx = mfi.tilebox();
      auto const& state_arr = ldata_p->state.array(mfi, FIRSTSOOT);
      SootData* sd = soot_model->getSootData_d();
      amrex::ParallelFor(
        gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          GpuArray<Real, NUM_SOOT_MOMENTS + 1> moments;
          for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; mom++) {
            moments[mom] = state_arr(i, j, k, mom);
          }
          sd->momConvClipConv(moments.data());
          for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; mom++) {
            state_arr(i, j, k, mom) = moments[mom];
          }
        });
    }
  }
}

// TODO: This isn't working yet
#if 0
void
PeleLM::addSootDerivePlotVars(PeleLMDeriveList& derive_lst)
{
  // Add in soot variables
  Vector<std::string> sootNames = {"rho_soot", "sum_rho_soot"};
  derive_lst.add(
    "soot_vars", IndexType::TheCellType(), sootNames.size(), sootNames,
    soot_genvars, PeleLMDeriveRec::TheSameBox);

  // Variables associated with the second mode (large particles)
  Vector<std::string> large_part_names = {"NL", "soot_V_L", "soot_S_L"};
  derive_lst.add(
    "soot_large_particles", IndexType::TheCellType(), large_part_names.size(),
    large_part_names, soot_largeparticledata, PeleLMDeriveRec::TheSameBox);
}
#endif
#endif
