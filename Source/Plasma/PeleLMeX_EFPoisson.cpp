#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <PeleLMeX_EF_Constants.H>
#include <PeleLMeX_DiffusionOp.H>

void
PeleLM::poissonSolveEF(const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::poissonSolveEF()");
  if (ef_verbose) {
    amrex::Print() << " EF Poisson solve \n";
  }

  // Get the phiV BCRec
  auto bcRecPhiV = fetchBCRecArray(PHIV, 1);

  // Build Poisson RHS: charge distribution
  int nGhost = 0;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> rhsPoisson(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    rhsPoisson[lev].reset(new amrex::MultiFab(
      grids[lev], dmap[lev], 1, nGhost, amrex::MFInfo(), *m_factory[lev]));

    auto ldata_p = getLevelDataPtr(lev, a_time);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*rhsPoisson[lev], amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& rhoY = ldata_p->state.const_array(mfi, FIRSTSPEC);
      auto const& nE = ldata_p->state.const_array(mfi, NE);
      auto const& rhs = rhsPoisson[lev]->array(mfi);
      amrex::Real factor = -1.0; // / ( eps0  * epsr);
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
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
    0, {}, {}, {}, 0, bcRecPhiV, 1, 1, -eps0 * epsr, {});
}
