#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <PeleLMeX_EF_Constants.H>
#include <AMReX_FillPatchUtil.H>
#include <PeleLMeX_BCfill.H>

using namespace amrex;

void
PeleLM::ionDriftVelocity(std::unique_ptr<AdvanceAdvData>& advData)
{
  //----------------------------------------------------------------
  // set udrift boundaries to zero
  if (advData->uDrift[0][0].nGrow() > 0) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      AMREX_D_TERM(advData->uDrift[lev][0].setBndry(0.0);
                   , advData->uDrift[lev][1].setBndry(0.0);
                   , advData->uDrift[lev][2].setBndry(0.0););
    }
  }

  //----------------------------------------------------------------
  // Get the gradient of Old and New phiV
  Vector<Array<MultiFab, AMREX_SPACEDIM>> gphiVOld(finest_level + 1);
  Vector<Array<MultiFab, AMREX_SPACEDIM>> gphiVNew(finest_level + 1);
  int nGrow = 0; // No need for ghost face on gphiV
  for (int lev = 0; lev <= finest_level; ++lev) {
    const auto& ba = grids[lev];
    const auto& factory = Factory(lev);
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      gphiVOld[lev][idim].define(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
        nGrow, MFInfo(), factory);
      gphiVNew[lev][idim].define(
        amrex::convert(ba, IntVect::TheDimensionVector(idim)), dmap[lev], 1,
        nGrow, MFInfo(), factory);
    }
  }

  int do_avgDown = 0; // TODO or should I ?
  auto bcRecPhiV = fetchBCRecArray(PHIV, 1);
  getDiffusionOp()->computeGradient(
    GetVecOfArrOfPtrs(gphiVOld), {}, // don't need the laplacian out
    GetVecOfConstPtrs(getPhiVVect(AmrOldTime)), bcRecPhiV[0], do_avgDown);
  getDiffusionOp()->computeGradient(
    GetVecOfArrOfPtrs(gphiVNew), {}, // don't need the laplacian out
    GetVecOfConstPtrs(getPhiVVect(AmrNewTime)), bcRecPhiV[0], do_avgDown);

  //----------------------------------------------------------------
  // TODO : this assumes that all the ions are grouped together at th end ...
  auto bcRecIons =
    fetchBCRecArray(FIRSTSPEC + NUM_SPECIES - NUM_IONS, NUM_IONS);

  for (int lev = 0; lev <= finest_level; ++lev) {
    // Get CC ions t^{n+1/2} mobilities
    // TODO In the old version, there's a switch to only use an instant. value
    auto ldataOld_p = getLevelDataPtr(lev, AmrOldTime);
    auto ldataNew_p = getLevelDataPtr(lev, AmrNewTime);

    MultiFab mobH_cc(grids[lev], dmap[lev], NUM_IONS, 1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mobH_cc, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& gbx = mfi.growntilebox();
      const auto& mob_o = ldataOld_p->mob_cc.const_array(mfi);
      const auto& mob_n = ldataNew_p->mob_cc.const_array(mfi);
      const auto& mob_h = mobH_cc.array(mfi);
      amrex::ParallelFor(
        gbx, NUM_IONS,
        [mob_o, mob_n,
         mob_h] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          mob_h(i, j, k, n) = 0.5 * (mob_o(i, j, k, n) + mob_n(i, j, k, n));
        });
    }

    // Get the face centered ions mobility
    int doZeroVisc = 0;
    Array<MultiFab, AMREX_SPACEDIM> mobH_ec =
      getDiffusivity(lev, 0, NUM_IONS, doZeroVisc, bcRecIons, mobH_cc);

    // Assemble the ions drift velocity
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(mobH_ec[idim], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.tilebox();
        const auto& mob_h = mobH_ec[idim].const_array(mfi);
        const auto& gp_o = gphiVOld[lev][idim].const_array(mfi);
        const auto& gp_n = gphiVNew[lev][idim].const_array(mfi);
        const auto& Ud_Sp = advData->uDrift[lev][idim].array(mfi);
        amrex::ParallelFor(
          bx, NUM_IONS,
          [mob_h, gp_o, gp_n,
           Ud_Sp] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            Ud_Sp(i, j, k, n) =
              mob_h(i, j, k, n) * -0.5 * (gp_o(i, j, k) + gp_n(i, j, k));
          });
      }
    }
  }

  //----------------------------------------------------------------
  // Average down faces
  for (int lev = finest_level; lev > 0; --lev) {
#ifdef AMREX_USE_EB
    EB_average_down_faces(
      GetArrOfConstPtrs(advData->uDrift[lev]),
      GetArrOfPtrs(advData->uDrift[lev - 1]), refRatio(lev - 1), geom[lev - 1]);
#else
    average_down_faces(
      GetArrOfConstPtrs(advData->uDrift[lev]),
      GetArrOfPtrs(advData->uDrift[lev - 1]), refRatio(lev - 1), geom[lev - 1]);
#endif
  }

  // FillPatch Udrift on levels > 0
  for (int lev = 0; lev <= finest_level; ++lev) {
    if (lev > 0) {
      IntVect rr = geom[lev].Domain().size() / geom[lev - 1].Domain().size();
      Interpolater* mapper = &face_linear_interp;

      // Set BCRec for Umac
      Vector<BCRec> bcrec(NUM_IONS);
      for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
        for (int ion = 0; ion < NUM_IONS; ion++) {
          if (geom[lev - 1].isPeriodic(idim)) {
            bcrec[ion].setLo(idim, BCType::int_dir);
            bcrec[ion].setHi(idim, BCType::int_dir);
          } else {
            bcrec[ion].setLo(idim, BCType::foextrap);
            bcrec[ion].setHi(idim, BCType::foextrap);
          }
        }
      }
      Array<Vector<BCRec>, AMREX_SPACEDIM> bcrecArr = {
        AMREX_D_DECL(bcrec, bcrec, bcrec)};

      PhysBCFunct<GpuBndryFuncFab<umacFill>> crse_bndry_func(
        geom[lev - 1], bcrec, umacFill{});
      Array<PhysBCFunct<GpuBndryFuncFab<umacFill>>, AMREX_SPACEDIM>
        cbndyFuncArr = {
          AMREX_D_DECL(crse_bndry_func, crse_bndry_func, crse_bndry_func)};

      PhysBCFunct<GpuBndryFuncFab<umacFill>> fine_bndry_func(
        geom[lev], bcrec, umacFill{});
      Array<PhysBCFunct<GpuBndryFuncFab<umacFill>>, AMREX_SPACEDIM>
        fbndyFuncArr = {
          AMREX_D_DECL(fine_bndry_func, fine_bndry_func, fine_bndry_func)};

      Real dummy = 0.;
      FillPatchTwoLevels(
        GetArrOfPtrs(advData->uDrift[lev]), IntVect(1), dummy,
        {GetArrOfPtrs(advData->uDrift[lev - 1])}, {dummy},
        {GetArrOfPtrs(advData->uDrift[lev])}, {dummy}, 0, 0, NUM_IONS,
        geom[lev - 1], geom[lev - 1], cbndyFuncArr, 0, fbndyFuncArr, 0, rr,
        mapper, bcrecArr, 0);
    } else {
      AMREX_D_TERM(
        advData->uDrift[lev][0].FillBoundary(geom[lev].periodicity());
        , advData->uDrift[lev][1].FillBoundary(geom[lev].periodicity());
        , advData->uDrift[lev][2].FillBoundary(geom[lev].periodicity()));
    }
  }
}

void
PeleLM::ionDriftAddUmac(int lev, std::unique_ptr<AdvanceAdvData>& advData)
{
  // Add umac to the ions drift velocity to get the effective velocity
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(advData->umac[lev][idim], TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const Box gbx = mfi.growntilebox();
      const auto& umac = advData->umac[lev][idim].const_array(mfi);
      const auto& Ud_Sp = advData->uDrift[lev][idim].array(mfi);
      amrex::ParallelFor(
        gbx, NUM_IONS,
        [umac, Ud_Sp] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          Ud_Sp(i, j, k, n) += umac(i, j, k);
        });
    }
    advData->uDrift[lev][idim].FillBoundary(geom[lev].periodicity());
  }
}
