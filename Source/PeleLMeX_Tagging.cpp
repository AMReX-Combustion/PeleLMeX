#include <PeleLMeX.H>
#include <AMReX_TagBox.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBAmrUtil.H>
#endif

void
PeleLM::ErrorEst(
  int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ng*/)
{
  BL_PROFILE("PeleLMeX::ErrorEst()");

#ifdef AMREX_USE_EB
  // Tag EB up to m_EB_refine_LevMax-1 if Static or
  //              m_EB_refine_LevAdapt-1 if Adaptive
  if (
    (m_EB_refine_type == "Static" && lev < m_EB_refine_LevMax) ||
    (m_EB_refine_type == "Adaptive" && lev < m_EB_refine_LevAdapt)) {
    const amrex::MultiFab& state = (getLevelDataPtr(lev, AmrNewTime))->state;
    TagCutCells(tags, state);
  }
#endif

  for (const auto& errTag : errTags) {
    std::unique_ptr<amrex::MultiFab> mf;
    if (!errTag.Field().empty()) {
      mf = deriveComp(errTag.Field(), time, lev, errTag.NGrow());
    }
    errTag(
      tags, mf.get(), amrex::TagBox::CLEAR, amrex::TagBox::SET, time, lev,
      geom[lev]);
  }

#ifdef AMREX_USE_EB
  // Untag covered cells
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(tags, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const auto& bx = mfi.tilebox();
    auto tag = tags.array(mfi);
    auto vfrac = EBFactory(lev).getVolFrac().const_array(mfi);
    amrex::ParallelFor(
      bx, [tag, vfrac] AMREX_GPU_HOST_DEVICE(int i, int j, int k) {
        if (vfrac(i, j, k) <= 0.0) {
          tag(i, j, k) = amrex::TagBox::CLEAR;
        }
      });
  }

  // Untag cell close to EB
  if (m_EB_refine_type == "Static" && lev >= m_EB_refine_LevMax) {
    // Get distance function at current level
    amrex::MultiFab signDist(
      grids[lev], dmap[lev], 1, 0, amrex::MFInfo(), EBFactory(lev));
    getEBDistance(lev, signDist);

    // Estimate how far I need to derefine
    constexpr amrex::Real sqrt2 = 1.4142135623730951;
    const amrex::Real diagFac = sqrt2 * m_derefineEBBuffer;
    amrex::Real clearTagDist =
      Geom(m_EB_refine_LevMax).CellSize(0) *
      static_cast<amrex::Real>(nErrorBuf(m_EB_refine_LevMax)) * diagFac;
    for (int ilev = m_EB_refine_LevMax + 1; ilev <= finest_level; ++ilev) {
      clearTagDist += static_cast<amrex::Real>(nErrorBuf(ilev)) *
                      Geom(m_EB_refine_LevMax).CellSize(0) * diagFac;
    }

    // Untag cells too close to EB
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(tags, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const auto& bx = mfi.tilebox();
      const auto& dist = signDist.const_array(mfi);
      auto tag = tags.array(mfi);
      amrex::ParallelFor(
        bx,
        [dist, clearTagDist, tag] AMREX_GPU_HOST_DEVICE(int i, int j, int k) {
          if (dist(i, j, k) < clearTagDist) {
            tag(i, j, k) = amrex::TagBox::CLEAR;
          }
        });
    }
  }
#endif
}
