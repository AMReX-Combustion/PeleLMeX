#include <AMReX_buildInfo.H>
#include <PeleLMeX.H>
#include <PeleLMeX_K.H>
#include <hydro_utils.H>
#include <memory>
#ifdef PELE_USE_EFIELD
#include <PeleLMeX_EF_Constants.H>
#endif

using namespace amrex;

void
writeBuildInfo()
{
  std::string OtherLine = std::string(78, '-') + "\n";
  std::string SkipSpace = std::string(8, ' ');

  // build information
  std::cout << PrettyLine;
  std::cout << " PeleLMeX Build Information\n";
  std::cout << PrettyLine;

  std::cout << "build date:    " << buildInfoGetBuildDate() << "\n";
  std::cout << "build machine: " << buildInfoGetBuildMachine() << "\n";
  std::cout << "build dir:     " << buildInfoGetBuildDir() << "\n";
  std::cout << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

  std::cout << "\n";

  std::cout << "COMP:          " << buildInfoGetComp() << "\n";
  std::cout << "COMP version:  " << buildInfoGetCompVersion() << "\n";

  std::cout << "\n";

  std::cout << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
  std::cout << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

  std::cout << "\n";

  std::cout << "Link flags:    " << buildInfoGetLinkFlags() << "\n";
  std::cout << "Libraries:     " << buildInfoGetLibraries() << "\n";

  std::cout << "\n";

  for (int n = 1; n <= buildInfoGetNumModules(); n++) {
    std::cout << buildInfoGetModuleName(n) << ": " << buildInfoGetModuleVal(n)
              << "\n";
  }

  std::cout << "\n";

  const char* githash1 = buildInfoGetGitHash(1);
  const char* githash2 = buildInfoGetGitHash(2);
  const char* githash3 = buildInfoGetGitHash(3);

  if (strlen(githash1) > 0) {
    std::cout << "PeleLMeX     git describe: " << githash1 << "\n";
  }
  if (strlen(githash2) > 0) {
    std::cout << "AMReX        git describe: " << githash2 << "\n";
  }
  if (strlen(githash3) > 0) {
    std::cout << "PelePhysics  git describe: " << githash3 << "\n";
  }

  const char* buildgithash = buildInfoGetBuildGitHash();
  const char* buildgitname = buildInfoGetBuildGitName();
  if (strlen(buildgithash) > 0) {
    std::cout << buildgitname << " git describe: " << buildgithash << "\n";
  }

  std::cout << "\n\n";
}

void
PeleLM::fluxDivergence(
  const Vector<MultiFab*>& a_divergence,
  int div_comp,
  const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& a_fluxes,
  int flux_comp,
  int ncomp,
  int intensiveFluxes,
  Real scale)
{
  BL_PROFILE("PeleLMeX::fluxDivergence()");
  if (intensiveFluxes != 0) { // Fluxes are intensive -> need area scaling in
                              // div
    for (int lev = 0; lev <= finest_level; ++lev) {
      intFluxDivergenceLevel(
        lev, *a_divergence[lev], div_comp, a_fluxes[lev], flux_comp, ncomp,
        scale);
    }
  } else { // Fluxes are extensive
    for (int lev = 0; lev <= finest_level; ++lev) {
      extFluxDivergenceLevel(
        lev, *a_divergence[lev], div_comp, a_fluxes[lev], flux_comp, ncomp,
        scale);
    }
  }
}

void
PeleLM::fluxDivergence(
  const Vector<MultiFab*>& a_divergence,
  int div_comp,
  const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& a_fluxes,
  int flux_comp,
  const Vector<MultiFab*>& a_EBfluxes,
  int ebflux_comp,
  int ncomp,
  int intensiveFluxes,
  Real scale)
{

  BL_PROFILE("PeleLMeX::fluxDivergence()");
  if (intensiveFluxes != 0) { // Fluxes are intensive -> need area scaling in
                              // div
    for (int lev = 0; lev <= finest_level; ++lev) {
      intFluxDivergenceLevelEB(
        lev, *a_divergence[lev], div_comp, a_fluxes[lev], flux_comp,
        a_EBfluxes[lev], ebflux_comp, ncomp, scale);
    }
  } else { // Fluxes are extensive
    for (int lev = 0; lev <= finest_level; ++lev) {
      extFluxDivergenceLevel(
        lev, *a_divergence[lev], div_comp, a_fluxes[lev], flux_comp, ncomp,
        scale);
    }
  }
}

void
PeleLM::fluxDivergenceRD(
  const Vector<const MultiFab*>& a_state,
  int state_comp,
  const Vector<MultiFab*>& a_divergence,
  int div_comp,
  const Vector<Array<MultiFab*, AMREX_SPACEDIM>>& a_fluxes,
  int flux_comp,
  const Vector<MultiFab*>& a_EBfluxes,
  int ebflux_comp,
  int ncomp,
  int intensiveFluxes,
  const BCRec* state_bc_d,
  const Real& scale,
  const Real& a_dt)
{
  BL_PROFILE("PeleLMeX::fluxDivergenceRD()");
#ifdef AMREX_USE_EB
  int have_ebfluxes = (a_EBfluxes.empty()) ? 0 : 1;
  for (int lev = 0; lev <= finest_level; ++lev) {
    //----------------------------------------------------------------
    // Use a temporary MF to hold divergence before redistribution
    int nGrow_divTmp = 3;
    MultiFab divTmp(
      grids[lev], dmap[lev], ncomp, nGrow_divTmp, MFInfo(), EBFactory(lev));
    divTmp.setVal(0.0);
    if (intensiveFluxes != 0) { // Fluxes are intensive -> need area scaling in
                                // div
      if (have_ebfluxes != 0) {
        intFluxDivergenceLevelEB(
          lev, divTmp, 0, a_fluxes[lev], flux_comp, a_EBfluxes[lev],
          ebflux_comp, ncomp, scale);
      } else {
        intFluxDivergenceLevel(
          lev, divTmp, 0, a_fluxes[lev], flux_comp, ncomp, scale);
      }
    } else { // Fluxes are extensive
      extFluxDivergenceLevel(
        lev, divTmp, 0, a_fluxes[lev], flux_comp, ncomp, scale);
    }

    // Need FillBoundary before redistribution
    divTmp.FillBoundary(geom[lev].periodicity());

    // Redistribute diffusion term
    redistributeDiff(
      lev, a_dt, divTmp, 0, *a_divergence[lev], div_comp, *a_state[lev],
      state_comp, ncomp, state_bc_d, geom[lev]);
  }
#else
  amrex::ignore_unused(a_state);
  amrex::ignore_unused(state_comp);
  amrex::ignore_unused(state_bc_d);
  amrex::ignore_unused(a_dt);
  amrex::ignore_unused(a_EBfluxes);
  amrex::ignore_unused(ebflux_comp);
  fluxDivergence(
    a_divergence, div_comp, a_fluxes, flux_comp, ncomp, intensiveFluxes, scale);
#endif
}

void
PeleLM::extFluxDivergenceLevel(
  int lev,
  MultiFab& a_divergence,
  int div_comp,
  const Array<MultiFab*, AMREX_SPACEDIM>& a_fluxes,
  int flux_comp,
  int ncomp,
  Real scale)
{

  AMREX_ASSERT(a_divergence.nComp() >= div_comp + ncomp);

  // Get the volume
  MultiFab volume(grids[lev], dmap[lev], 1, 0);
  geom[lev].GetVolume(volume);

#ifdef AMREX_USE_EB
  auto const& ebfact = EBFactory(lev);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(a_divergence, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    AMREX_D_TERM(auto const& fluxX = a_fluxes[0]->const_array(mfi, flux_comp);
                 , auto const& fluxY = a_fluxes[1]->const_array(mfi, flux_comp);
                 ,
                 auto const& fluxZ = a_fluxes[2]->const_array(mfi, flux_comp););
    auto const& divergence = a_divergence.array(mfi, div_comp);
    auto const& vol = volume.const_array(mfi);

#ifdef AMREX_USE_EB
    auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    auto const& flag = flagfab.const_array();
#endif

#ifdef AMREX_USE_EB
    if (flagfab.getType(bx) == FabType::covered) { // Covered boxes
      amrex::ParallelFor(
        bx, ncomp,
        [divergence] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          divergence(i, j, k, n) = 0.0;
        });
    } else if (flagfab.getType(bx) != FabType::regular) { // EB containing boxes
      auto vfrac = ebfact.getVolFrac().const_array(mfi);
      amrex::ParallelFor(
        bx, [ncomp, flag, vfrac, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ),
             vol, scale] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (flag(i, j, k).isCovered()) {
            for (int n = 0; n < ncomp; n++) {
              divergence(i, j, k, n) = 0.0;
            }
          } else if (flag(i, j, k).isRegular()) {
            extFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale,
              divergence);
          } else {
            Real vfracinv = 1.0 / vfrac(i, j, k);
            extFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale,
              divergence);
            for (int n = 0; n < ncomp; n++) {
              divergence(i, j, k, n) *= vfracinv;
            }
          }
        });
    } else // Regular boxes
#endif
    {
      amrex::ParallelFor(
        bx, [ncomp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol,
             scale] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          extFluxDivergence_K(
            i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ), vol, scale,
            divergence);
        });
    }
  }
}

void
PeleLM::intFluxDivergenceLevel(
  int lev,
  MultiFab& a_divergence,
  int div_comp,
  const Array<MultiFab*, AMREX_SPACEDIM>& a_fluxes,
  int flux_comp,
  int ncomp,
  Real scale)
{

  AMREX_ASSERT(a_divergence.nComp() >= div_comp + ncomp);

  // Get the volume
  MultiFab volume(grids[lev], dmap[lev], 1, 0);
  geom[lev].GetVolume(volume);

  // Get areas
  const Real* dx = Geom(lev).CellSize();
#if (AMREX_SPACEDIM == 2)
  MultiFab mf_ax, mf_ay;
  if (geom[lev].IsRZ()) {
    geom[lev].GetFaceArea(mf_ax, grids[lev], dmap[lev], 0, 0);
    geom[lev].GetFaceArea(mf_ay, grids[lev], dmap[lev], 1, 0);
  }
  Real areax = dx[1];
  Real areay = dx[0];
#elif (AMREX_SPACEDIM == 3)
  Real areax = dx[1] * dx[2];
  Real areay = dx[0] * dx[2];
  Real areaz = dx[0] * dx[1];
#endif

  // Get areafrac if EB
#ifdef AMREX_USE_EB
  auto const& ebfact = EBFactory(lev);
  Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
  areafrac = ebfact.getAreaFrac();
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(a_divergence, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    AMREX_D_TERM(auto const& fluxX = a_fluxes[0]->const_array(mfi, flux_comp);
                 , auto const& fluxY = a_fluxes[1]->const_array(mfi, flux_comp);
                 ,
                 auto const& fluxZ = a_fluxes[2]->const_array(mfi, flux_comp););
    auto const& divergence = a_divergence.array(mfi, div_comp);
    auto const& vol = volume.const_array(mfi);

#ifdef AMREX_USE_EB
    auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    auto const& flag = flagfab.const_array();

    if (flagfab.getType(bx) == FabType::covered) { // Covered boxes
      amrex::ParallelFor(
        bx, ncomp,
        [divergence] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          divergence(i, j, k, n) = 0.0;
        });
    } else if (flagfab.getType(bx) != FabType::regular) { // EB containing boxes
      auto vfrac = ebfact.getVolFrac().const_array(mfi);
      AMREX_D_TERM(const auto& afrac_x = areafrac[0]->array(mfi);
                   , const auto& afrac_y = areafrac[1]->array(mfi);
                   , const auto& afrac_z = areafrac[2]->array(mfi););
      amrex::ParallelFor(
        bx, [ncomp, flag, vfrac, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ),
             AMREX_D_DECL(afrac_x, afrac_y, afrac_z),
             AMREX_D_DECL(areax, areay, areaz), vol,
             scale] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (flag(i, j, k).isCovered()) {
            for (int n = 0; n < ncomp; n++) {
              divergence(i, j, k, n) = 0.0;
            }
          } else if (flag(i, j, k).isRegular()) {
            intFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ),
              AMREX_D_DECL(areax, areay, areaz), vol, scale, divergence);
          } else {
            Real vfracinv = 1.0 / vfrac(i, j, k);
            EB_intFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ),
              AMREX_D_DECL(afrac_x, afrac_y, afrac_z),
              AMREX_D_DECL(areax, areay, areaz), vol, scale, divergence);
            for (int n = 0; n < ncomp; n++) {
              divergence(i, j, k, n) *= vfracinv;
            }
          }
        });
    } else // Regular boxes
#endif
    {
#if (AMREX_SPACEDIM == 2)
      if (geom[lev].IsRZ()) {
        Array4<Real const> const& ax = mf_ax.const_array(mfi);
        Array4<Real const> const& ay = mf_ay.const_array(mfi);
        amrex::ParallelFor(
          bx, [ncomp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ), ax, ay,
               vol, scale] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            intFluxDivergence_rz_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ), ax, ay, vol,
              scale, divergence);
          });
      } else
#endif
      {
        amrex::ParallelFor(
          bx, [ncomp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ),
               AMREX_D_DECL(areax, areay, areaz), vol,
               scale] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            intFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ),
              AMREX_D_DECL(areax, areay, areaz), vol, scale, divergence);
          });
      }
    }
  }
}

void
PeleLM::intFluxDivergenceLevelEB(
  int lev,
  MultiFab& a_divergence,
  int div_comp,
  const Array<MultiFab*, AMREX_SPACEDIM>& a_fluxes,
  int flux_comp,
  const MultiFab* a_EBfluxes,
  int ebflux_comp,
  int ncomp,
  Real scale)
{

  AMREX_ASSERT(a_divergence.nComp() >= div_comp + ncomp);

  // Get the volume
  MultiFab volume(grids[lev], dmap[lev], 1, 0);
  geom[lev].GetVolume(volume);

  // Get area
  const GpuArray<Real, AMREX_SPACEDIM> dx = Geom(lev).CellSizeArray();
#if (AMREX_SPACEDIM == 2)
  Real areax = dx[1];
  Real areay = dx[0];
#elif (AMREX_SPACEDIM == 3)
  Real areax = dx[1] * dx[2];
  Real areay = dx[0] * dx[2];
  Real areaz = dx[0] * dx[1];
#endif

  // Get areafrac if EB
#ifdef AMREX_USE_EB
  auto const& ebfact = EBFactory(lev);
  Array<const MultiCutFab*, AMREX_SPACEDIM> areafrac;
  areafrac = ebfact.getAreaFrac();
  const auto* eb_area = &(ebfact.getBndryArea());
#else
  amrex::ignore_unused(a_EBfluxes, ebflux_comp);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(a_divergence, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.tilebox();
    AMREX_D_TERM(auto const& fluxX = a_fluxes[0]->const_array(mfi, flux_comp);
                 , auto const& fluxY = a_fluxes[1]->const_array(mfi, flux_comp);
                 ,
                 auto const& fluxZ = a_fluxes[2]->const_array(mfi, flux_comp););
    auto const& divergence = a_divergence.array(mfi, div_comp);
    auto const& vol = volume.const_array(mfi);

#ifdef AMREX_USE_EB
    auto const& ebflux = a_EBfluxes->const_array(mfi, ebflux_comp);
    auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    auto const& flag = flagfab.const_array();

    if (flagfab.getType(bx) == FabType::covered) { // Covered boxes
      amrex::ParallelFor(
        bx, ncomp,
        [divergence] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          divergence(i, j, k, n) = 0.0;
        });
    } else if (flagfab.getType(bx) != FabType::regular) { // EB containing boxes
      auto vfrac = ebfact.getVolFrac().const_array(mfi);
      AMREX_D_TERM(const auto& afrac_x = areafrac[0]->array(mfi);
                   , const auto& afrac_y = areafrac[1]->array(mfi);
                   , const auto& afrac_z = areafrac[2]->array(mfi););
      const auto& ebarea = eb_area->array(mfi);
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          if (flag(i, j, k).isCovered()) {
            for (int n = 0; n < ncomp; n++) {
              divergence(i, j, k, n) = 0.0;
            }
          } else if (flag(i, j, k).isRegular()) {
            intFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ),
              AMREX_D_DECL(areax, areay, areaz), vol, scale, divergence);
          } else {
            Real vfracinv = 1.0 / vfrac(i, j, k);
            EB_intFluxDivergence_K(
              i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ),
              AMREX_D_DECL(afrac_x, afrac_y, afrac_z),
              AMREX_D_DECL(areax, areay, areaz), ebflux, ebarea, vol, dx[0],
              scale, divergence);
            for (int n = 0; n < ncomp; n++) {
              divergence(i, j, k, n) *= vfracinv;
            }
          }
        });
    } else // Regular boxes
#endif
    {
      amrex::ParallelFor(
        bx, [ncomp, divergence, AMREX_D_DECL(fluxX, fluxY, fluxZ),
             AMREX_D_DECL(areax, areay, areaz), vol,
             scale] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          intFluxDivergence_K(
            i, j, k, ncomp, AMREX_D_DECL(fluxX, fluxY, fluxZ),
            AMREX_D_DECL(areax, areay, areaz), vol, scale, divergence);
        });
    }
  }
}

void
PeleLM::
  advFluxDivergence( // NOLINT(readability-convert-member-functions-to-static)
    int a_lev,
    MultiFab& a_divergence,
    int div_comp,
    MultiFab& a_divu,
    const Array<const MultiFab*, AMREX_SPACEDIM>& a_fluxes,
    int flux_comp,
    const Array<const MultiFab*, AMREX_SPACEDIM>& a_faceState,
    int face_comp,
    int ncomp,
    int const* l_conserv_d,
    const Geometry& a_geom,
    amrex::Real scale,
    bool fluxes_are_area_weighted) const
{
  BL_PROFILE("PeleLMeX::advFluxDivergence()");

  AMREX_ASSERT(a_divergence.nComp() >= div_comp + ncomp);
  AMREX_ASSERT(a_fluxes[0]->nComp() >= flux_comp + ncomp);
  AMREX_ASSERT(a_faceState[0]->nComp() >= face_comp + ncomp);

#if (AMREX_SPACEDIM == 2)
  MultiFab mf_ax, mf_ay;
  if (geom[a_lev].IsRZ()) {
    geom[a_lev].GetFaceArea(mf_ax, grids[a_lev], dmap[a_lev], 0, 0);
    geom[a_lev].GetFaceArea(mf_ay, grids[a_lev], dmap[a_lev], 1, 0);
  }
#else
  amrex::ignore_unused(a_lev);
#endif

#ifdef AMREX_USE_EB
  auto const& ebfact = EBFactory(a_lev);
#endif

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(a_divergence, TilingIfNotGPU()); mfi.isValid(); ++mfi) {

    Box const& bx = mfi.tilebox();

    // Get the divergence
    auto const& div_arr = a_divergence.array(mfi, div_comp);
    AMREX_D_TERM(auto const& fx = a_fluxes[0]->const_array(mfi, flux_comp);
                 , auto const& fy = a_fluxes[1]->const_array(mfi, flux_comp);
                 , auto const& fz = a_fluxes[2]->const_array(mfi, flux_comp);)

#ifdef AMREX_USE_EB
    auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
    auto const& vfrac_arr = ebfact.getVolFrac().const_array(mfi);
    if (flagfab.getType(bx) == FabType::singlevalued) {
      HydroUtils::EB_ComputeDivergence(
        bx, div_arr, AMREX_D_DECL(fx, fy, fz), vfrac_arr, ncomp, a_geom, scale,
        fluxes_are_area_weighted);
    } else if (flagfab.getType(bx) == FabType::regular)
#endif
    {
      HydroUtils::ComputeDivergence(
        bx, div_arr, AMREX_D_DECL(fx, fy, fz), ncomp, a_geom, scale,
        fluxes_are_area_weighted);
    }

    // If convective, we define u dot grad q = div (u q) - q div(u)
    // averaging face and t^{n+1/2} q to the cell center
    auto const& divu_arr = a_divu.const_array(mfi);
    AMREX_D_TERM(
      auto const& facex = a_faceState[0]->const_array(mfi, face_comp);
      , auto const& facey = a_faceState[1]->const_array(mfi, face_comp);
      , auto const& facez = a_faceState[2]->const_array(mfi, face_comp);)

#ifdef AMREX_USE_EB
    if (flagfab.getType(bx) == FabType::covered) {
      AMREX_PARALLEL_FOR_4D(
        bx, ncomp, i, j, k, n, { div_arr(i, j, k, n) = 0.0; });
    } else if (flagfab.getType(bx) == FabType::singlevalued) {
      AMREX_D_TERM(
        auto const& apx_arr = ebfact.getAreaFrac()[0]->const_array(mfi);
        , auto const& apy_arr = ebfact.getAreaFrac()[1]->const_array(mfi);
        , auto const& apz_arr = ebfact.getAreaFrac()[2]->const_array(mfi););
      ParallelFor(
        bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
          if ((l_conserv_d[n] == 0) && vfrac_arr(i, j, k) > 0.) {
            Real qwsum = AMREX_D_TERM(
              apx_arr(i, j, k) * facex(i, j, k, n) +
                apx_arr(i + 1, j, k) * facex(i + 1, j, k, n),
              +apy_arr(i, j, k) * facey(i, j, k, n) +
                apy_arr(i, j + 1, k) * facey(i, j + 1, k, n),
              +apz_arr(i, j, k) * facez(i, j, k, n) +
                apz_arr(i, j, k + 1) * facez(i, j, k + 1, n));
            Real areasum = AMREX_D_TERM(
              apx_arr(i, j, k) + apx_arr(i + 1, j, k),
              +apy_arr(i, j, k) + apy_arr(i, j + 1, k),
              +apz_arr(i, j, k) + apz_arr(i, j, k + 1));
            // Note that because we define adv update as MINUS div(u q), here we
            // add q div (u)
            div_arr(i, j, k, n) += qwsum / areasum * divu_arr(i, j, k);
          }
        });
    } else
#endif
    {
#if (AMREX_SPACEDIM == 2)
      if (geom[a_lev].IsRZ()) {
        Array4<Real const> const& ax = mf_ax.const_array(mfi);
        Array4<Real const> const& ay = mf_ay.const_array(mfi);
        ParallelFor(
          bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (l_conserv_d[n] == 0) {
              Real qavg = AMREX_D_TERM(
                ax(i, j, k) * facex(i, j, k, n) +
                  ax(i + 1, j, k) * facex(i + 1, j, k, n),
                +ay(i, j, k) * facey(i, j, k, n) +
                  ay(i, j + 1, k) * facey(i, j + 1, k, n),
                +0.0);
              Real areasum =
                ax(i, j, k) + ax(i + 1, j, k) + ay(i, j, k) + ay(i, j + 1, k);
              qavg /= areasum;
              // Note that because we define adv update as MINUS div(u q), here
              // we add q div (u)
              div_arr(i, j, k, n) += qavg * divu_arr(i, j, k);
            }
          });
      } else
#endif
      {
        ParallelFor(
          bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            if (l_conserv_d[n] == 0) {
              Real qavg = AMREX_D_TERM(
                facex(i, j, k, n) + facex(i + 1, j, k, n),
                +facey(i, j, k, n) + facey(i, j + 1, k, n),
                +facez(i, j, k, n) + facez(i, j, k + 1, n));
              AMREX_D_PICK(qavg *= 0.5;, qavg *= 0.25;, qavg /= 6.0;)
              // Note that because we define adv update as MINUS div(u q), here
              // we add q div (u)
              div_arr(i, j, k, n) += qavg * divu_arr(i, j, k);
            }
          });
      }
    }
  }
}

void
PeleLM::floorSpecies(const TimeStamp& a_time)
{
  BL_PROFILE("PeleLMeX::floorSpecies()");
  AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);
  if (m_floor_species == 0) {
    return;
  }

  for (int lev = 0; lev <= finest_level; ++lev) {

    auto* ldata_p = getLevelDataPtr(lev, a_time);
    auto const& sma = ldata_p->state.arrays();

    amrex::ParallelFor(
      ldata_p->state,
      [=] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
        fabMinMax(
          i, j, k, NUM_SPECIES, 0.0, AMREX_REAL_MAX,
          Array4<Real>(sma[box_no], FIRSTSPEC));
#ifdef PELE_USE_EFIELD
        fabMinMax(
          i, j, k, 1, 0.0, AMREX_REAL_MAX, Array4<Real>(sma[box_no], NE));
#endif
        // Update density accordingly ...
        sma[box_no](i, j, k, DENSITY) = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
          sma[box_no](i, j, k, DENSITY) += sma[box_no](i, j, k, FIRSTSPEC + n);
        }

        // ... as well as rhoh
        auto eos = pele::physics::PhysicsType::eos();
        Real massfrac[NUM_SPECIES] = {0.0};
        Real rhoinv = Real(1.0) / sma[box_no](i, j, k, DENSITY);
        for (int n = 0; n < NUM_SPECIES; n++) {
          massfrac[n] = sma[box_no](i, j, k, FIRSTSPEC + n) * rhoinv;
        }
        Real h_cgs = 0.0;
        eos.TY2H(sma[box_no](i, j, k, TEMP), massfrac, h_cgs);
        sma[box_no](i, j, k, RHOH) =
          h_cgs * 1.0e-4 * sma[box_no](i, j, k, DENSITY);
      });
    Gpu::streamSynchronize();
  }
}

void
PeleLM::resetCoveredMask()
{
  BL_PROFILE("PeleLMeX::resetCoveredMask()");
  if (m_resetCoveredMask != 0) {

    if (m_verbose != 0) {
      Print() << " Resetting fine-covered cells mask \n";
    }

    for (int lev = 0; lev < finest_level; ++lev) {
      // Set a fine-covered mask
      BoxArray baf = grids[lev + 1];
      baf.coarsen(ref_ratio[lev]);
      m_coveredMask[lev]->setVal(1);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      {
        std::vector<std::pair<int, Box>> isects;
        for (MFIter mfi(*m_coveredMask[lev], TilingIfNotGPU()); mfi.isValid();
             ++mfi) {
          auto const& mask = m_coveredMask[lev]->array(mfi);
          baf.intersections(grids[lev][mfi.index()], isects);
          for (const auto& is : isects) {
            amrex::ParallelFor(
              is.second, [mask] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                mask(i, j, k) = 0;
              });
          }
        }
      }

      //----------------------------------------------------------------------------
      // Setup a BoxArray for the chemistry

      // Get an uncovered BoxArray
      BoxArray baCompDom = complementIn(geom[lev].Domain(), baf);
      BoxArray baUnCovered = intersect(baCompDom, grids[lev]);
      // Chop in smaller boxes if triggered
      if (m_max_grid_size_chem.min() > 0) {
        baUnCovered.maxSize(m_max_grid_size_chem);
      }

      // Assemble a BoxArray with covered and uncovered ones + flags
      BoxList bl(grids[lev].ixType());
      bl.reserve(baUnCovered.size() + baf.size());
      m_baChemFlag[lev].resize(baUnCovered.size() + baf.size());
      int bxIdx = 0;
      for (int i = 0, n = static_cast<int>(baUnCovered.size()); i < n;
           ++i) { // Append uncovered boxes
        bl.push_back(baUnCovered[i]);
        m_baChemFlag[lev][bxIdx] = 1;
        bxIdx += 1;
      }
      for (int i = 0, n = static_cast<int>(baf.size()); i < n;
           ++i) { // Append covered boxes
        bl.push_back(baf[i]);
        m_baChemFlag[lev][bxIdx] = 0;
        bxIdx += 1;
      }
      m_baChem[lev] = std::make_unique<BoxArray>(std::move(bl));
      m_dmapChem[lev] = std::make_unique<DistributionMapping>(*m_baChem[lev]);

      // Load balancing of the chemistry DMap
      if (m_doLoadBalance != 0) {
        loadBalanceChemLev(lev);
      }
    }

    // Set a BoxArray for the chemistry on the finest level too
    m_baChem[finest_level] = std::make_unique<BoxArray>(grids[finest_level]);
    if (m_max_grid_size_chem.min() > 0) {
      m_baChem[finest_level]->maxSize(m_max_grid_size_chem);
    }
    m_baChemFlag[finest_level].resize(m_baChem[finest_level]->size());
    std::fill(
      m_baChemFlag[finest_level].begin(), m_baChemFlag[finest_level].end(), 1);
    m_dmapChem[finest_level] =
      std::make_unique<DistributionMapping>(*m_baChem[finest_level]);

    if ((m_doLoadBalance != 0) && m_max_grid_size_chem.min() > 0) {
      loadBalanceChemLev(finest_level);
    }

    // Switch off trigger
    m_resetCoveredMask = 0;

  } else {
    // Just load balance the chem. distribution map
    if (m_doLoadBalance != 0) {
      loadBalanceChem();
    }
  }

  //----------------------------------------------------------------------------
  // Need to compute the uncovered volume
  if (m_uncoveredVol < 0.0) {
    Vector<MultiFab> dummy(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
      dummy[lev].define(grids[lev], dmap[lev], 1, 0, MFInfo(), *m_factory[lev]);
      dummy[lev].setVal(1.0);
    }
    m_uncoveredVol = MFSum(GetVecOfConstPtrs(dummy), 0);
  }
}

void
PeleLM::loadBalanceChem()
{

  for (int lev = 0; lev <= finest_level; ++lev) {
    // Finest grid uses AmrCore DM unless different max grid size specified.
    // Keep the AmrCore DM.
    if (lev == finest_level && m_max_grid_size_chem.min() < 0) {
      continue;
    }
    loadBalanceChemLev(lev);
  }
}

void
PeleLM::loadBalanceChemLev(int a_lev)
{

  LayoutData<Real> new_cost(*m_baChem[a_lev], *m_dmapChem[a_lev]);
  computeCosts(a_lev, new_cost, m_loadBalanceCostChem);

  // Use efficiency: average MPI rank cost / max cost
  amrex::Real currentEfficiency = 0.0;
  amrex::Real testEfficiency = 0.0;

  DistributionMapping test_dmap;
  // Build the test dmap, w/o braodcasting
  if (m_loadBalanceMethodChem == LoadBalanceMethod::SFC) {

    test_dmap = DistributionMapping::makeSFC(
      new_cost, currentEfficiency, testEfficiency, false,
      ParallelDescriptor::IOProcessorNumber());

  } else if (m_loadBalanceMethodChem == LoadBalanceMethod::Knapsack) {

    const amrex::Real navg = static_cast<Real>(m_baChem[a_lev]->size()) /
                             static_cast<Real>(ParallelDescriptor::NProcs());
    const int nmax = static_cast<int>(
      std::max(std::round(m_loadBalanceKSfactor * navg), std::ceil(navg)));
    test_dmap = DistributionMapping::makeKnapSack(
      new_cost, currentEfficiency, testEfficiency, nmax, false,
      ParallelDescriptor::IOProcessorNumber());
  }

  // IO proc determine if the test dmap offers significant improvements
  int updateDmap = 0;
  if (
    (m_loadBalanceEffRatioThreshold > 0.0) &&
    (ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber())) {
    updateDmap = static_cast<int>(
      testEfficiency > m_loadBalanceEffRatioThreshold * currentEfficiency);
  }
  ParallelDescriptor::Bcast(
    &updateDmap, 1, ParallelDescriptor::IOProcessorNumber());

  if (m_verbose > 2 && (updateDmap != 0)) {
    Print() << "   Old Chem LoadBalancing efficiency on a_lev " << a_lev << ": "
            << currentEfficiency << "\n"
            << "   New Chem LoadBalancing efficiency: " << testEfficiency
            << " \n";
  }

  // Bcast the test dmap if better
  if (updateDmap != 0) {
    Vector<int> pmap;
    if (
      ParallelDescriptor::MyProc() == ParallelDescriptor::IOProcessorNumber()) {
      pmap = test_dmap.ProcessorMap();
    } else {
      pmap.resize(static_cast<std::size_t>(m_baChem[a_lev]->size()));
    }
    ParallelDescriptor::Bcast(
      pmap.data(), pmap.size(), ParallelDescriptor::IOProcessorNumber());

    if (
      ParallelDescriptor::MyProc() != ParallelDescriptor::IOProcessorNumber()) {
      test_dmap = DistributionMapping(pmap);
    }
    m_dmapChem[a_lev] = std::make_unique<DistributionMapping>(test_dmap);
  }
}

// Return a unique_ptr with the entire derive
std::unique_ptr<MultiFab>
PeleLM::derive(const std::string& a_name, Real a_time, int lev, int nGrow)
{
  BL_PROFILE("PeleLMeX::derive()");
  AMREX_ASSERT(nGrow >= 0);

  std::unique_ptr<MultiFab> mf;

  bool itexists = derive_lst.canDerive(a_name) || isStateVariable(a_name) ||
                  isReactVariable(a_name);

  if (!itexists) {
    amrex::Error("PeleLM::derive(): unknown variable: " + a_name);
  }

  const PeleLMDeriveRec* rec = derive_lst.get(a_name);

  if (rec != nullptr) { // This is a derived variable
    mf = std::make_unique<MultiFab>(
      grids[lev], dmap[lev], rec->numDerive(), nGrow, MFInfo(), Factory(lev));
    std::unique_ptr<MultiFab> statemf =
      fillPatchState(lev, a_time, m_nGrowState);
    // Get pressure: TODO no fillpatch for pressure just yet, simply get new
    // state
    auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
    std::unique_ptr<MultiFab> reactmf;
    if (m_do_react != 0) {
      reactmf = fillPatchReact(lev, a_time, nGrow);
    }
    auto stateBCs = fetchBCRecArray(VELX, NVAR);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(nGrow);
      FArrayBox& derfab = (*mf)[mfi];
      FArrayBox const& statefab = (*statemf)[mfi];
      FArrayBox const& reactfab =
        (m_incompressible) != 0 ? ldata_p->press[mfi] : (*reactmf)[mfi];
      FArrayBox const& pressfab = ldata_p->press[mfi];
      rec->derFunc()(
        this, bx, derfab, 0, rec->numDerive(), statefab, reactfab, pressfab,
        geom[lev], a_time, stateBCs, lev);
    }
  } else if (isStateVariable(a_name)) { // This is a state variable
    mf = std::make_unique<MultiFab>(
      grids[lev], dmap[lev], 1, nGrow, MFInfo(), Factory(lev));
    int idx = stateVariableIndex(a_name);
    std::unique_ptr<MultiFab> statemf = fillPatchState(lev, a_time, nGrow);
    MultiFab::Copy(*mf, *statemf, idx, 0, 1, nGrow);
  } else { // This is a reaction variable
    mf = std::make_unique<MultiFab>(
      grids[lev], dmap[lev], 1, nGrow, MFInfo(), Factory(lev));
    int idx = reactVariableIndex(a_name);
    std::unique_ptr<MultiFab> reactmf = fillPatchReact(lev, a_time, nGrow);
    MultiFab::Copy(*mf, *reactmf, idx, 0, 1, nGrow);
  }

  return mf;
}

// Return a unique_ptr with only the required component of a derive
std::unique_ptr<MultiFab>
PeleLM::deriveComp(const std::string& a_name, Real a_time, int lev, int nGrow)
{
  BL_PROFILE("PeleLMeX::derive()");
  AMREX_ASSERT(nGrow >= 0);

  std::unique_ptr<MultiFab> mf;

  bool itexists = derive_lst.canDerive(a_name) || isStateVariable(a_name) ||
                  isReactVariable(a_name);

  if (!itexists) {
    amrex::Error("PeleLM::derive(): unknown variable: " + a_name);
  }

  const PeleLMDeriveRec* rec = derive_lst.get(a_name);

  if (rec != nullptr) { // This is a derived variable
    mf = std::make_unique<MultiFab>(
      grids[lev], dmap[lev], 1, nGrow, MFInfo(), Factory(lev));
    std::unique_ptr<MultiFab> statemf =
      fillPatchState(lev, a_time, m_nGrowState);
    // Get pressure: TODO no fillpatch for pressure just yet, simply get new
    // state
    auto* ldata_p = getLevelDataPtr(lev, AmrNewTime);
    std::unique_ptr<MultiFab> reactmf = fillPatchReact(lev, a_time, nGrow);
    auto stateBCs = fetchBCRecArray(VELX, NVAR);

    // Temp MF for all the derive components
    MultiFab derTemp(grids[lev], dmap[lev], rec->numDerive(), nGrow);
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*mf, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.growntilebox(nGrow);
      FArrayBox& derfab = derTemp[mfi];
      FArrayBox const& statefab = (*statemf)[mfi];
      FArrayBox const& reactfab =
        (m_incompressible) != 0 ? ldata_p->press[mfi] : (*reactmf)[mfi];
      FArrayBox const& pressfab = ldata_p->press[mfi];
      rec->derFunc()(
        this, bx, derfab, 0, rec->numDerive(), statefab, reactfab, pressfab,
        geom[lev], a_time, stateBCs, lev);
    }
    // Copy into outgoing unique_ptr
    int derComp = rec->variableComp(a_name);
    if (derComp < 0) {
      amrex::Error(
        "PeleLM::deriveComp(): unknown derive component: " + a_name + " of " +
        rec->variableName(1000));
    }
    MultiFab::Copy(*mf, derTemp, derComp, 0, 1, nGrow);
  } else if (isStateVariable(a_name)) { // This is a state variable
    mf = std::make_unique<MultiFab>(
      grids[lev], dmap[lev], 1, nGrow, MFInfo(), Factory(lev));
    int idx = stateVariableIndex(a_name);
    std::unique_ptr<MultiFab> statemf = fillPatchState(lev, a_time, nGrow);
    MultiFab::Copy(*mf, *statemf, idx, 0, 1, nGrow);
  } else { // This is a reaction variable
    mf = std::make_unique<MultiFab>(
      grids[lev], dmap[lev], 1, nGrow, MFInfo(), Factory(lev));
    int idx = reactVariableIndex(a_name);
    std::unique_ptr<MultiFab> reactmf = fillPatchReact(lev, a_time, nGrow);
    MultiFab::Copy(*mf, *reactmf, idx, 0, 1, nGrow);
  }

  return mf;
}

void
PeleLM::initProgressVariable()
{
  Vector<std::string> varNames;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    varNames);
  varNames.push_back("temp");

  ParmParse pp("peleLM");
  std::string Cformat;
  int hasUserC = static_cast<int>(pp.contains("progressVariable.format"));
  if (hasUserC != 0) {
    pp.query("progressVariable.format", Cformat);
    if (Cformat == "Cantera") { // use a Cantera-like format with
                                // <entry>:<weight>, default to 0.0
      // Weights
      Vector<std::string> stringIn;
      Vector<Real> weightsIn(NUM_SPECIES + 1, 0.0);
      int entryCount = pp.countval("progressVariable.weights");
      stringIn.resize(entryCount);
      pp.getarr("progressVariable.weights", stringIn, 0, entryCount);
      parseVars(varNames, stringIn, weightsIn);

      // Cold side/Hot side
      Vector<Real> coldState(NUM_SPECIES + 1, 0.0);
      entryCount = pp.countval("progressVariable.coldState");
      stringIn.resize(entryCount);
      pp.getarr("progressVariable.coldState", stringIn, 0, entryCount);
      parseVars(varNames, stringIn, coldState);
      Vector<Real> hotState(NUM_SPECIES + 1, 0.0);
      entryCount = pp.countval("progressVariable.hotState");
      stringIn.resize(entryCount);
      pp.getarr("progressVariable.hotState", stringIn, 0, entryCount);
      parseVars(varNames, stringIn, hotState);
      m_C0 = 0.0;
      m_C1 = 0.0;
      for (int i = 0; i < NUM_SPECIES + 1; ++i) {
        m_Cweights[i] = weightsIn[i];
        m_C0 += coldState[i] * m_Cweights[i];
        m_C1 += hotState[i] * m_Cweights[i];
      }
    } else if (Cformat == "RealList") { // use a list of Real. MUST
                                        // contains an entry
                                        // for each species+Temp
      // Weights
      Vector<Real> weightsIn;
      int entryCount = pp.countval("progressVariable.weights");
      AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES + 1);
      weightsIn.resize(entryCount);
      pp.getarr("progressVariable.weights", weightsIn, 0, entryCount);
      for (int i = 0; i < NUM_SPECIES; ++i) {
        m_Cweights[i] = weightsIn[i];
      }
      // Cold side/Hot side
      entryCount = pp.countval("progressVariable.coldState");
      AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES + 1);
      Vector<Real> coldState(entryCount);
      pp.getarr("progressVariable.coldState", coldState, 0, entryCount);
      entryCount = pp.countval("progressVariable.hotState");
      AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES + 1);
      Vector<Real> hotState(entryCount);
      pp.getarr("progressVariable.hotState", hotState, 0, entryCount);
      m_C0 = 0.0;
      m_C1 = 0.0;
      for (int i = 0; i < NUM_SPECIES + 1; ++i) {
        m_C0 += coldState[i] * m_Cweights[i];
        m_C1 += hotState[i] * m_Cweights[i];
      }
    } else {
      Abort(
        "Unknown progressVariable.format ! Should be 'Cantera' or 'RealList'");
    }
    pp.query("progressVariable.revert", m_Crevert);
  }
}

void
PeleLM::parseVars(
  const Vector<std::string>& a_varsNames,
  const Vector<std::string>& a_stringIn,
  Vector<Real>& a_rVars)
{
  const int varCountIn = static_cast<int>(a_stringIn.size());

  // For each entry in the user-provided composition, parse name and value
  std::string delimiter = ":";
  for (int i = 0; i < varCountIn; i++) {
    long unsigned sep = a_stringIn[i].find(delimiter);
    if (sep == std::string::npos) {
      Abort(
        "Error parsing '" + a_stringIn[i] + "' --> unable to find delimiter :");
    }
    std::string varNameIn = a_stringIn[i].substr(0, sep);
    Real value =
      std::stod(a_stringIn[i].substr(sep + 1, a_stringIn[i].length()));
    int foundIt = 0;
    for (int k = 0; k < a_varsNames.size(); k++) {
      if (varNameIn == a_varsNames[k]) {
        a_rVars[k] = value;
        foundIt = 1;
      }
    }
    if (foundIt == 0) {
      Abort(
        "Error parsing '" + a_stringIn[i] +
        "' --> unable to match to any provided variable name");
    }
  }
}

Real
PeleLM::MLNorm0(const Vector<const MultiFab*>& a_MF)
{
  BL_PROFILE("PeleLMeX::MLNorm0()");
  Real r = 0.0;
  for (int lev = 0; lev < a_MF.size(); ++lev) {
    if (lev != finest_level) {
      r = std::max(r, a_MF[lev]->norm0(*m_coveredMask[lev], 0, 0, true));
    } else {
      r = std::max(r, a_MF[lev]->norm0(0, 0, true, true));
    }
  }
  ParallelDescriptor::ReduceRealMax(r);
  return r;
}

Vector<Real>
PeleLM::MLNorm0(const Vector<const MultiFab*>& a_MF, int startcomp, int ncomp)
{
  BL_PROFILE("PeleLMeX::MLNorm0()");
  AMREX_ASSERT(a_MF[0]->nComp() >= startcomp + ncomp);
  Vector<Real> r(ncomp);
  for (int n = 0; n < ncomp; n++) {
    r[n] = 0.0;
  }
  for (int lev = 0; lev < a_MF.size(); ++lev) {
    if (lev != finest_level) {
      for (int n = 0; n < ncomp; n++) {
        r[n] = std::max(
          r[n], a_MF[lev]->norm0(*m_coveredMask[lev], startcomp + n, 0, true));
      }
    } else {
      for (int n = 0; n < ncomp; n++) {
        r[n] = std::max(r[n], a_MF[lev]->norm0(startcomp + n, 0, true, true));
      }
    }
  }
  ParallelDescriptor::ReduceRealMax(r.data(), ncomp);
  return r;
}

bool
PeleLM::isStateVariable(std::string_view a_name)
{
  // Check state
  return std::any_of(
    stateComponents.begin(), stateComponents.end(),
    [=](const auto& stateComponent) {
      return std::get<1>(stateComponent) == a_name;
    });
}

bool
PeleLM::isReactVariable(std::string_view a_name)
{
  // Check reaction state
  return std::any_of(
    reactComponents.begin(), reactComponents.end(),
    [=](const auto& reactComponent) {
      return std::get<1>(reactComponent) == a_name;
    });
}

int
PeleLM::stateVariableIndex(std::string_view a_name)
{
  int idx = -1;
  if (!isStateVariable(a_name)) {
    amrex::Error(
      "PeleLM::stateVariableIndex(): unknown State variable: " +
      static_cast<std::string>(a_name));
  }
  for (const auto& stateComponent : stateComponents) {
    if (std::get<1>(stateComponent) == a_name) {
      idx = std::get<0>(stateComponent);
    }
  }
  return idx;
}

int
PeleLM::reactVariableIndex(std::string_view a_name)
{
  int idx = -1;
  if (!isReactVariable(a_name)) {
    amrex::Error(
      "PeleLM::reactVariableIndex(): unknown Reaction variable: " +
      static_cast<std::string>(a_name));
  }
  for (const auto& reactComponent : reactComponents) {
    if (std::get<1>(reactComponent) == a_name) {
      idx = std::get<0>(reactComponent);
    }
  }
  return idx;
}

Vector<int>
PeleLM::fetchAdvTypeArray(int scomp, int ncomp)
{
  Vector<int> types(ncomp);
  for (int comp = 0; comp < ncomp; comp++) {
    types[comp] = m_AdvTypeState[scomp + comp];
  }
  return types;
}

Vector<int>
PeleLM::fetchDiffTypeArray(int scomp, int ncomp)
{
  Vector<int> types(ncomp);
  for (int comp = 0; comp < ncomp; comp++) {
    types[comp] = m_DiffTypeState[scomp + comp];
  }
  return types;
}

Real
PeleLM::MFSum(const Vector<const MultiFab*>& a_mf, int comp)
{
  BL_PROFILE("PeleLMeX::MFSum()");
  // Get the integral of the MF, not including the fine-covered and
  // EB-covered cells

  Real volwgtsum = 0.0;

  for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef AMREX_USE_EB
    // For EB, use constant vol
    const Real* dx = geom[lev].CellSize();
    Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    // Use amrex::ReduceSum
    auto const& ebfact = dynamic_cast<EBFArrayBoxFactory const&>(Factory(lev));
    auto const& vfrac = ebfact.getVolFrac();

    Real sm = 0.0;
    if (lev != finest_level) {
      sm = amrex::ReduceSum(
        *a_mf[lev], vfrac, *m_coveredMask[lev], 0,
        [vol, comp] AMREX_GPU_HOST_DEVICE(
          Box const& bx, Array4<Real const> const& mf_arr,
          Array4<Real const> const& vf_arr,
          Array4<int const> const& covered_arr) -> Real {
          Real sum = 0.0;
          AMREX_LOOP_3D(bx, i, j, k, {
            sum += mf_arr(i, j, k, comp) * vf_arr(i, j, k) * vol *
                   static_cast<Real>(covered_arr(i, j, k));
          });
          return sum;
        });
    } else {
      sm = amrex::ReduceSum(
        *a_mf[lev], vfrac, 0,
        [vol, comp] AMREX_GPU_HOST_DEVICE(
          Box const& bx, Array4<Real const> const& mf_arr,
          Array4<Real const> const& vf_arr) -> Real {
          Real sum = 0.0;
          AMREX_LOOP_3D(bx, i, j, k, {
            sum += mf_arr(i, j, k, comp) * vf_arr(i, j, k) * vol;
          });
          return sum;
        });
    }
#else
    // Get the geometry volume to account for 2D-RZ
    MultiFab volume(grids[lev], dmap[lev], 1, 0);
    geom[lev].GetVolume(volume);

    Real sm = 0.0;
    if (lev != finest_level) {
      sm = amrex::ReduceSum(
        *a_mf[lev], volume, *m_coveredMask[lev], 0,
        [comp] AMREX_GPU_HOST_DEVICE(
          Box const& bx, Array4<Real const> const& mf_arr,
          Array4<Real const> const& vol_arr,
          Array4<int const> const& covered_arr) -> Real {
          Real sum = 0.0;
          AMREX_LOOP_3D(bx, i, j, k, {
            sum += mf_arr(i, j, k, comp) * vol_arr(i, j, k) *
                   static_cast<Real>(covered_arr(i, j, k));
          });
          return sum;
        });
    } else {
      sm = amrex::ReduceSum(
        *a_mf[lev], volume, 0,
        [comp] AMREX_GPU_HOST_DEVICE(
          Box const& bx, Array4<Real const> const& mf_arr,
          Array4<Real const> const& vol_arr) -> Real {
          Real sum = 0.0;
          AMREX_LOOP_3D(
            bx, i, j, k, { sum += mf_arr(i, j, k, comp) * vol_arr(i, j, k); });
          return sum;
        });
    }
#endif

    volwgtsum += sm;
  } // lev

  ParallelDescriptor::ReduceRealSum(volwgtsum);

  return volwgtsum;
}

/*
Array<Real,3>
PeleLM::MFStat (const Vector<const MultiFab*> &a_mf, int comp)
{
   // Get the min/max/mean of a given component, not including the fine-covered
cells
}
*/

void
PeleLM::setTypicalValues(const TimeStamp& a_time, int is_init)
{
  // Get state Max/Min
  auto stateMax =
    (m_incompressible) != 0
      ? MLmax(GetVecOfConstPtrs(getStateVect(a_time)), 0, AMREX_SPACEDIM)
      : MLmax(GetVecOfConstPtrs(getStateVect(a_time)), 0, NVAR);
  auto stateMin =
    (m_incompressible) != 0
      ? MLmin(GetVecOfConstPtrs(getStateVect(a_time)), 0, AMREX_SPACEDIM)
      : MLmin(GetVecOfConstPtrs(getStateVect(a_time)), 0, NVAR);

  // Fill typical values vector
  for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
    typical_values[idim] =
      std::max(stateMax[VELX + idim], std::abs(stateMin[VELX + idim]));
  }

  if (m_incompressible == 0) {
    // Average between max/min
    typical_values[DENSITY] = 0.5 * (stateMax[DENSITY] + stateMin[DENSITY]);
    for (int n = 0; n < NUM_SPECIES; n++) {
      typical_values[FIRSTSPEC + n] =
        0.5 * (stateMax[FIRSTSPEC + n] + stateMin[FIRSTSPEC + n]) /
        typical_values[DENSITY];
    }
    typical_values[RHOH] =
      0.5 * (stateMax[RHOH] + stateMin[RHOH]) / typical_values[DENSITY];
    typical_values[TEMP] = 0.5 * (stateMax[TEMP] + stateMin[TEMP]);
    typical_values[RHORT] = m_pOld;
#ifdef PELE_USE_EFIELD
    typical_values[NE] = 0.5 * (stateMax[NE] + stateMin[NE]);
#endif

    // Pass into chemsitry if requested
    updateTypicalValuesChem();
  }

  if ((is_init != 0) || m_verbose > 1) {
    Print() << PrettyLine;
    Print() << " Typical values: " << '\n';
    Print() << "\tVelocity: ";
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      Print() << typical_values[idim] << ' ';
    }
    Print() << '\n';
    if (m_incompressible == 0) {
      Print() << "\tDensity:  " << typical_values[DENSITY] << '\n';
      Print() << "\tTemp:     " << typical_values[TEMP] << '\n';
      Print() << "\tH:        " << typical_values[RHOH] << '\n';
      Vector<std::string> spec_names;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
        spec_names);
      for (int n = 0; n < NUM_SPECIES; n++) {
        Print() << "\tY_" << spec_names[n]
                << std::setw(
                     std::max(0, static_cast<int>(8 - spec_names[n].length())))
                << std::left << ":" << typical_values[FIRSTSPEC + n] << '\n';
      }
#ifdef PELE_USE_EFIELD
      Print() << "\tnE:       " << typical_values[NE] << '\n';
#endif
    }
    Print() << PrettyLine;
  }
}

void
PeleLM::updateTypicalValuesChem()
{
  if ((m_useTypValChem != 0) && (m_do_react != 0)) {
    if (m_verbose > 2) {
      Print() << " Update chemistry typical values \n";
    }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
      Vector<Real> typical_values_chem;
      typical_values_chem.resize(NUM_SPECIES + 1);
      for (int i = 0; i < NUM_SPECIES; ++i) {
        typical_values_chem[i] = amrex::max(
          m_typicalYvalMin * typical_values[DENSITY] * 1.E-3,
          typical_values[FIRSTSPEC + i] * typical_values[DENSITY] *
            1.E-3); // CGS -> MKS conversion
      }
      typical_values_chem[NUM_SPECIES] = typical_values[TEMP];
#ifdef PELE_USE_EFIELD
      auto eos = pele::physics::PhysicsType::eos();
      Real mwt[NUM_SPECIES] = {0.0};
      eos.molecular_weight(mwt);
      typical_values_chem[E_ID] =
        typical_values[NE] / Na * mwt[E_ID] * 1.0e-6 * 1.0e-2;
#endif
      m_reactor->set_typ_vals_ode(typical_values_chem);
    }
  }
}

// MultiFab max, excluding EB-covered/fine-covered cells, local
Real
PeleLM::MFmax(const MultiFab* a_MF, const iMultiFab& a_mask, int comp)
{
  BL_PROFILE("PeleLMeX::MFmax()");
  Real mx = std::numeric_limits<Real>::lowest();

#ifdef AMREX_USE_EB
  if (a_MF->hasEBFabFactory()) {
    const auto& ebfactory =
      dynamic_cast<EBFArrayBoxFactory const&>(a_MF->Factory());
    auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
      auto const& flagsma = flags.const_arrays();
      auto const& ma = a_MF->const_arrays();
      auto const& mask = a_mask.const_arrays();
      mx = ParReduce(
        TypeList<ReduceOpMax>{}, TypeList<Real>{}, *a_MF, IntVect(0),
        [=] AMREX_GPU_DEVICE(
          int box_no, int i, int j, int k) noexcept -> GpuTuple<Real> {
          if (flagsma[box_no](i, j, k).isCovered() || !mask[box_no](i, j, k)) {
            return AMREX_REAL_LOWEST;
          } else {
            return ma[box_no](i, j, k, comp);
          }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max : mx)
#endif
      for (MFIter mfi(*a_MF, true); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        if (flags[mfi].getType(bx) != FabType::covered) {
          auto const& flag = flags.const_array(mfi);
          auto const& a = a_MF->const_array(mfi);
          auto const& mask = a_mask.const_array(mfi);
          AMREX_LOOP_3D(bx, i, j, k, {
            if (!flag(i, j, k).isCovered() && mask(i, j, k)) {
              mx = std::max(mx, a(i, j, k, comp));
            }
          });
        }
      }
    }
  } else
#endif
  {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
      auto const& ma = a_MF->const_arrays();
      auto const& mask = a_mask.const_arrays();
      mx = ParReduce(
        TypeList<ReduceOpMax>{}, TypeList<Real>{}, *a_MF, IntVect(0),
        [=] AMREX_GPU_DEVICE(
          int box_no, int i, int j, int k) noexcept -> GpuTuple<Real> {
          if (!mask[box_no](i, j, k)) {
            return AMREX_REAL_LOWEST;
          } else {
            return ma[box_no](i, j, k, comp);
          }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max : mx)
#endif
      for (MFIter mfi(*a_MF, true); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        auto const& a = a_MF->const_array(mfi);
        auto const& mask = a_mask.const_array(mfi);
        AMREX_LOOP_3D(bx, i, j, k, {
          if (mask(i, j, k)) {
            mx = std::max(mx, a(i, j, k, comp));
          }
        });
      }
    }
  }

  return mx;
}

// MultiFab min, excluding EB-covered/fine-covered cells, local
Real
PeleLM::MFmin(const MultiFab* a_MF, const iMultiFab& a_mask, int comp)
{
  BL_PROFILE("PeleLMeX::MFmin()");
  Real mn = std::numeric_limits<Real>::max();

#ifdef AMREX_USE_EB
  if (a_MF->hasEBFabFactory()) {
    const auto& ebfactory =
      dynamic_cast<EBFArrayBoxFactory const&>(a_MF->Factory());
    auto const& flags = ebfactory.getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
      auto const& flagsma = flags.const_arrays();
      auto const& ma = a_MF->const_arrays();
      auto const& mask = a_mask.const_arrays();
      mn = ParReduce(
        TypeList<ReduceOpMin>{}, TypeList<Real>{}, *a_MF, IntVect(0),
        [=] AMREX_GPU_DEVICE(
          int box_no, int i, int j, int k) noexcept -> GpuTuple<Real> {
          if (flagsma[box_no](i, j, k).isCovered() || !mask[box_no](i, j, k)) {
            return AMREX_REAL_MAX;
          } else {
            return ma[box_no](i, j, k, comp);
          }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min : mn)
#endif
      for (MFIter mfi(*a_MF, true); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        if (flags[mfi].getType(bx) != FabType::covered) {
          auto const& flag = flags.const_array(mfi);
          auto const& a = a_MF->const_array(mfi);
          auto const& mask = a_mask.const_array(mfi);
          AMREX_LOOP_3D(bx, i, j, k, {
            if (!flag(i, j, k).isCovered() && mask(i, j, k)) {
              mn = std::min(mn, a(i, j, k, comp));
            }
          });
        }
      }
    }
  } else
#endif
  {
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
      auto const& ma = a_MF->const_arrays();
      auto const& mask = a_mask.const_arrays();
      mn = ParReduce(
        TypeList<ReduceOpMin>{}, TypeList<Real>{}, *a_MF, IntVect(0),
        [=] AMREX_GPU_DEVICE(
          int box_no, int i, int j, int k) noexcept -> GpuTuple<Real> {
          if (!mask[box_no](i, j, k)) {
            return AMREX_REAL_MAX;
          } else {
            return ma[box_no](i, j, k, comp);
          }
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min : mn)
#endif
      for (MFIter mfi(*a_MF, true); mfi.isValid(); ++mfi) {
        Box const& bx = mfi.tilebox();
        auto const& a = a_MF->const_array(mfi);
        auto const& mask = a_mask.const_array(mfi);
        AMREX_LOOP_3D(bx, i, j, k, {
          if (mask(i, j, k)) {
            mn = std::min(mn, a(i, j, k, comp));
          }
        });
      }
    }
  }

  return mn;
}

// MultiLevel max, exlucing EB-covered/fine-covered cells
Vector<Real>
PeleLM::MLmax(const Vector<const MultiFab*>& a_MF, int scomp, int ncomp)
{
  BL_PROFILE("PeleLMeX::MLmax()");
  AMREX_ASSERT(a_MF[0]->nComp() >= scomp + ncomp);

  Vector<Real> nmax(ncomp, AMREX_REAL_LOWEST);

  for (int lev = 0; lev < a_MF.size(); ++lev) {
    if (lev != finest_level) {
      for (int n = 0; n < ncomp; n++) {
        nmax[n] =
          std::max(nmax[n], MFmax(a_MF[lev], *m_coveredMask[lev], scomp + n));
      }
    } else {
      for (int n = 0; n < ncomp; n++) {
        nmax[n] = std::max(nmax[n], a_MF[lev]->max(scomp + n, 0, true));
      }
    }
  }

  ParallelDescriptor::ReduceRealMax(nmax.data(), ncomp);
  return nmax;
}

// MultiLevel min, exlucing EB-covered/fine-covered cells
Vector<Real>
PeleLM::MLmin(const Vector<const MultiFab*>& a_MF, int scomp, int ncomp)
{
  BL_PROFILE("PeleLMeX::MLmin()");
  AMREX_ASSERT(a_MF[0]->nComp() >= scomp + ncomp);

  Vector<Real> nmin(ncomp, AMREX_REAL_MAX);

  for (int lev = 0; lev < a_MF.size(); ++lev) {
    if (lev != finest_level) {
      for (int n = 0; n < ncomp; n++) {
        nmin[n] =
          std::min(nmin[n], MFmin(a_MF[lev], *m_coveredMask[lev], scomp + n));
      }
    } else {
      for (int n = 0; n < ncomp; n++) {
        nmin[n] = std::min(nmin[n], a_MF[lev]->min(scomp + n, 0, true));
      }
    }
  }

  ParallelDescriptor::ReduceRealMin(nmin.data(), ncomp);
  return nmin;
}

void
PeleLM::checkMemory(const std::string& a_message) const
{
  if (m_checkMem == 0) {
    return;
  }

  const int IOProc = ParallelDescriptor::IOProcessorNumber();
#ifdef AMREX_USE_GPU
  Long free_mem_avail = Gpu::Device::freeMemAvailable() / (1024 * 1024);
  ParallelDescriptor::ReduceLongMin(free_mem_avail, IOProc);
  Print() << "     [" << a_message << "] GPU mem. avail. (MB) "
          << free_mem_avail << "\n";
#else
  // MultiFab memory usage
  Long max_fab_megabytes =
    amrex::TotalBytesAllocatedInFabsHWM() / (1024 * 1024);
  ParallelDescriptor::ReduceLongMax(max_fab_megabytes, IOProc);
  Print() << "     [" << a_message << "] MFs mem. allocated (MB) "
          << max_fab_megabytes << "\n";
#endif
}

void
PeleLM::initMixtureFraction()
{
  // Get default fuel and oxy tank composition: pure fuel vs air
  Vector<std::string> specNames;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    specNames);
  amrex::Real YF[NUM_SPECIES], YO[NUM_SPECIES];
  for (int i = 0; i < NUM_SPECIES; ++i) {
    YF[i] = 0.0;
    YO[i] = 0.0;
    if (specNames[i] == "O2") {
      YO[i] = 0.233;
    }
    if (specNames[i] == "N2") {
      YO[i] = 0.767;
    }
    if (i == fuelID) {
      YF[i] = 1.0;
    }
  }

  auto eos = pele::physics::PhysicsType::eos();
  // Overwrite with user-defined value if provided in input file
  ParmParse pp("peleLM");
  std::string MFformat;
  int hasUserMF = static_cast<int>(pp.contains("mixtureFraction.format"));
  if (hasUserMF != 0) {
    pp.query("mixtureFraction.format", MFformat);
    if (MFformat == "Cantera") { // use a Cantera-like format with
      // <SpeciesName>:<Value>, default in 0.0
      std::string MFCompoType;
      pp.query("mixtureFraction.type", MFCompoType);
      Vector<std::string> compositionIn;
      int entryCount = pp.countval("mixtureFraction.oxidTank");
      compositionIn.resize(entryCount);
      pp.getarr("mixtureFraction.oxidTank", compositionIn, 0, entryCount);
      parseComposition(compositionIn, MFCompoType, YO);
      entryCount = pp.countval("mixtureFraction.fuelTank");
      compositionIn.resize(entryCount);
      pp.getarr("mixtureFraction.fuelTank", compositionIn, 0, entryCount);
      parseComposition(compositionIn, MFCompoType, YF);
    } else if (MFformat == "RealList") { // use a list of Real. MUST
                                         // contains an entry
      // for each species in the mixture
      std::string MFCompoType;
      pp.query("mixtureFraction.type", MFCompoType);
      if (MFCompoType == "mass") {
        int entryCount = pp.countval("mixtureFraction.oxidTank");
        AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES);
        Vector<amrex::Real> compositionIn(NUM_SPECIES);
        pp.getarr("mixtureFraction.oxidTank", compositionIn, 0, NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; ++i) {
          YO[i] = compositionIn[i];
        }
        entryCount = pp.countval("mixtureFraction.fuelTank");
        AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES);
        pp.getarr("mixtureFraction.fuelTank", compositionIn, 0, NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; ++i) {
          YF[i] = compositionIn[i];
        }
      } else if (MFCompoType == "mole") {
        amrex::Real XF[NUM_SPECIES], XO[NUM_SPECIES];
        int entryCount = pp.countval("mixtureFraction.oxidTank");
        AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES);
        Vector<amrex::Real> compositionIn(NUM_SPECIES);
        pp.getarr("mixtureFraction.oxidTank", compositionIn, 0, NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; ++i) {
          XO[i] = compositionIn[i];
        }
        entryCount = pp.countval("mixtureFraction.fuelTank");
        AMREX_ALWAYS_ASSERT(entryCount == NUM_SPECIES);
        pp.getarr("mixtureFraction.fuelTank", compositionIn, 0, NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; ++i) {
          XF[i] = compositionIn[i];
        }

        eos.X2Y(XO, YO);
        eos.X2Y(XF, YF);
      } else {
        Abort("Unknown mixtureFraction.type ! Should be 'mass' or 'mole'");
      }
    } else {
      Abort(
        "Unknown mixtureFraction.format ! Should be 'Cantera' or 'RealList'");
    }
  }
  if (fuelID < 0 && (hasUserMF != 0)) {
    Print() << " Mixture fraction definition lacks fuelID: consider using "
               "peleLM.fuel_name keyword \n";
  }

  // Only interested in CHON -in that order. Compute Bilger weights
  amrex::Real atwCHON[4] = {0.0};
  pele::physics::eos::atomic_weightsCHON<pele::physics::PhysicsType::eos_type>(
    atwCHON);
  Beta_mix[0] = (atwCHON[0] != 0.0) ? 2.0 / atwCHON[0] : 0.0;
  Beta_mix[1] = (atwCHON[1] != 0.0) ? 1.0 / (2.0 * atwCHON[1]) : 0.0;
  Beta_mix[2] = (atwCHON[2] != 0.0) ? -1.0 / atwCHON[2] : 0.0;
  Beta_mix[3] = 0.0;

  // Compute each species weight for the Bilger formulation based on elemental
  // compo Only interested in CHON -in that order.
  int ecompCHON[NUM_SPECIES * 4];
  pele::physics::eos::element_compositionCHON<
    pele::physics::PhysicsType::eos_type>(ecompCHON);
  amrex::Real mwt[NUM_SPECIES];
  eos.molecular_weight(mwt);
  Zfu = 0.0;
  Zox = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i) {
    spec_Bilger_fact[i] = 0.0;
    for (int k = 0; k < 4; k++) {
      spec_Bilger_fact[i] +=
        Beta_mix[k] * (ecompCHON[i * 4 + k] * atwCHON[k] / mwt[i]);
    }
    Zfu += spec_Bilger_fact[i] * YF[i];
    Zox += spec_Bilger_fact[i] * YO[i];
  }
}

void
PeleLM::parseComposition(
  Vector<std::string> compositionIn,
  std::string compositionType,
  Real* massFrac)
{
  Real compoIn[NUM_SPECIES] = {0.0};

  // Get species names
  Vector<std::string> specNames;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    specNames);

  // For each entry in the user-provided composition, parse name and value
  std::string delimiter = ":";
  const int specCountIn = static_cast<int>(compositionIn.size());
  for (int i = 0; i < specCountIn; i++) {
    long unsigned sep = compositionIn[i].find(delimiter);
    if (sep == std::string::npos) {
      Abort(
        "Error parsing '" + compositionIn[i] +
        "' --> unable to find delimiter :");
    }
    std::string specNameIn = compositionIn[i].substr(0, sep);
    Real value =
      std::stod(compositionIn[i].substr(sep + 1, compositionIn[i].length()));
    int foundIt = 0;
    for (int k = 0; k < NUM_SPECIES; k++) {
      if (specNameIn == specNames[k]) {
        compoIn[k] = value;
        foundIt = 1;
      }
    }
    if (foundIt == 0) {
      Abort(
        "Error parsing '" + compositionIn[i] +
        "' --> unable to match to any species name");
    }
  }

  // Ensure that it sums to 1.0:
  Real sum = 0.0;
  for (double k : compoIn) {
    sum += k;
  }
  for (double& k : compoIn) {
    k /= sum;
  }

  // Fill the massFrac array, convert from mole fraction if necessary
  if (compositionType == "mass") { // mass
    for (int i = 0; i < NUM_SPECIES; i++) {
      massFrac[i] = compoIn[i];
    }
  } else if (compositionType == "mole") { // mole
    auto eos = pele::physics::PhysicsType::eos();
    eos.X2Y(compoIn, massFrac);
  } else {
    Abort("Unknown mixtureFraction.type ! Should be 'mass' or 'mole'");
  }
}

#ifdef AMREX_USE_EB
// Extend the cell-centered based signed distance function
void
PeleLM::extendSignedDistance(MultiFab* a_signDist, Real a_extendFactor)
{
  BL_PROFILE("PeleLMeX::extendSignedDistance()");
  // This is a not-so-pretty piece of code that'll take AMReX cell-averaged
  // signed distance and propagates it manually up to the point where we need to
  // have it for derefining.
  const auto geomdata = geom[0].data();
  Real maxSignedDist = a_signDist->max(0);
  const auto& ebfactory =
    dynamic_cast<EBFArrayBoxFactory const&>(a_signDist->Factory());
  const auto& flags = ebfactory.getMultiEBCellFlagFab();
  int nGrowFac = flags.nGrow() + 1;

  // First set the region far away at the max value we need
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
  for (MFIter mfi(*a_signDist, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const Box& bx = mfi.growntilebox();
    auto const& sd_cc = a_signDist->array(mfi);
    ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (sd_cc(i, j, k) >= maxSignedDist - 1e-12) {
        const Real* dx = geomdata.CellSize();
        sd_cc(i, j, k) = nGrowFac * dx[0] * a_extendFactor;
      }
    });
  }

  // Iteratively compute the distance function in boxes, propagating across
  // boxes using ghost cells If needed, increase the number of loop to extend
  // the reach of the distance function
  int nMaxLoop = 4;
  for (int dloop = 1; dloop <= nMaxLoop; dloop++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*a_signDist, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
      const Box& bx = mfi.tilebox();
      const Box& gbx = grow(bx, 1);
      if (flags[mfi].getType(gbx) == FabType::covered) {
        continue;
      }
      auto const& sd_cc = a_signDist->array(mfi);
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const auto glo = amrex::lbound(gbx);
        const auto ghi = amrex::ubound(gbx);
        const Real* dx = geomdata.CellSize();
        Real extendedDist = dx[0] * a_extendFactor;
        if (sd_cc(i, j, k) >= maxSignedDist - 1e-12) {
          Real closestEBDist = 1e12;
          for (int kk = glo.z; kk <= ghi.z; ++kk) {
            for (int jj = glo.y; jj <= ghi.y; ++jj) {
              for (int ii = glo.x; ii <= ghi.x; ++ii) {
                if ((i != ii) || (j != jj) || (k != kk)) {
                  if (sd_cc(ii, jj, kk) > 0.0) {
                    Real distToCell = std::sqrt(AMREX_D_TERM(
                      ((i - ii) * dx[0] * (i - ii) * dx[0]),
                      +((j - jj) * dx[1] * (j - jj) * dx[1]),
                      +((k - kk) * dx[2] * (k - kk) * dx[2])));
                    Real distToEB = distToCell + sd_cc(ii, jj, kk);
                    if (distToEB < closestEBDist) {
                      closestEBDist = distToEB;
                    }
                  }
                }
              }
            }
          }
          if (closestEBDist < 1e10) {
            sd_cc(i, j, k) = closestEBDist;
          } else {
            sd_cc(i, j, k) = extendedDist;
          }
        }
      });
    }
    a_signDist->FillBoundary(geom[0].periodicity());
  }
}
#endif
