#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <PeleLMeX_EF_Constants.H>
#include <PeleLMeX_DiffusionOp.H>

using namespace amrex;

void
PeleLM::poissonSolveEF(const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::poissonSolveEF()");
  if (ef_verbose) {
    Print() << " EF Poisson solve \n";
  }

  // Get the phiV BCRec
  auto bcRecPhiV = fetchBCRecArray(PHIV, 1);

  // Build Poisson RHS: charge distribution
  int nGhost = 0;
  Vector<std::unique_ptr<MultiFab>> rhsPoisson(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    rhsPoisson[lev].reset(new MultiFab(
      grids[lev], dmap[lev], 1, nGhost, MFInfo(), *m_factory[lev]));

    auto ldata_p = getLevelDataPtr(lev, a_time);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*rhsPoisson[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      auto const& rhoY = ldata_p->state.const_array(mfi, FIRSTSPEC);
      auto const& nE = ldata_p->state.const_array(mfi, NE);
      auto const& rhs = rhsPoisson[lev]->array(mfi);
      Real factor = -1.0; // / ( eps0  * epsr);
      amrex::ParallelFor(
        bx, [rhs, rhoY, nE, factor,
             zk = zk] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          rhs(i, j, k) = -nE(i, j, k) * elemCharge * factor;
          for (int n = 0; n < NUM_SPECIES; n++) {
            rhs(i, j, k) += zk[n] * rhoY(i, j, k, n) * factor;
          }
        });
    }
  }

  // Solve for PhiV
  getDiffusionOp()->diffuse_scalar(
    GetVecOfPtrs(getPhiVVect(a_time)), 0, GetVecOfConstPtrs(rhsPoisson), 0, {},
    0, {}, {}, {}, 0, bcRecPhiV, 1, 1, -eps0 * epsr);
}
