#include "PeleLMeX_Index.H"
#include "PeleLMeX.H"
#include "PeleLMeX_K.H"
#include "PeleLMeX_DeriveFunc.H"

#include <PelePhysics.H>
#include <mechanism.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBFabFactory.H>
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

//
// Extract temp
//
void
pelelmex_dertemp(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  amrex::ignore_unused(a_pelelm, ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_ASSERT(!a_pelelm->m_incompressible);
  auto const in_dat = statefab.array();
  auto der = derfab.array(dcomp);
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    der(i, j, k) = in_dat(i, j, k, TEMP);
  });
}

//
// Compute heat release
//
void
pelelmex_derheatrelease(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& reactfab,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  amrex::ignore_unused(a_pelelm, ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_ASSERT(!a_pelelm->m_incompressible);

  FArrayBox EnthFab;
  EnthFab.resize(bx, NUM_SPECIES, The_Async_Arena());

  auto const temp = statefab.const_array(TEMP);
  auto const react = reactfab.const_array(0);
  auto const& Hi = EnthFab.array();
  auto HRR = derfab.array(dcomp);
  auto const* leosparm = a_pelelm->eos_parms.device_parm();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    getHGivenT(i, j, k, temp, Hi, leosparm);
    HRR(i, j, k) = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      HRR(i, j, k) -= Hi(i, j, k, n) * react(i, j, k, n);
    }
  });
}

//
// Extract species mass fractions Y_n
//
void
pelelmex_dermassfrac(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  amrex::ignore_unused(a_pelelm, ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES + 1);
  AMREX_ASSERT(ncomp == NUM_SPECIES);
  AMREX_ASSERT(!a_pelelm->m_incompressible);
  auto const in_dat = statefab.array();
  auto der = derfab.array(dcomp);
  amrex::ParallelFor(
    bx, NUM_SPECIES, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      amrex::Real rhoinv = 1.0 / in_dat(i, j, k, DENSITY);
      der(i, j, k, n) = in_dat(i, j, k, FIRSTSPEC + n) * rhoinv;
    });
}

//
// Extract species mole fractions X_n
//
void
pelelmex_dermolefrac(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::ignore_unused(a_pelelm, ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES + 1);
  AMREX_ASSERT(ncomp == NUM_SPECIES);
  AMREX_ASSERT(!a_pelelm->m_incompressible);
  auto const in_dat = statefab.array();
  auto der = derfab.array(dcomp);
  auto const* leosparm = a_pelelm->eos_parms.device_parm();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real Yt[NUM_SPECIES] = {0.0};
    amrex::Real Xt[NUM_SPECIES] = {0.0};
    amrex::Real rhoinv = 1.0 / in_dat(i, j, k, DENSITY);
    for (int n = 0; n < NUM_SPECIES; n++) {
      Yt[n] = in_dat(i, j, k, FIRSTSPEC + n) * rhoinv;
    }
    auto eos = pele::physics::PhysicsType::eos(leosparm);
    eos.Y2X(Yt, Xt);
    for (int n = 0; n < NUM_SPECIES; n++) {
      der(i, j, k, n) = Xt[n];
    }
  });
}

//
// Extract rho - sum rhoY
//
void
pelelmex_derrhomrhoy(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  amrex::ignore_unused(a_pelelm, ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES + 1);
  AMREX_ASSERT(ncomp == 1);
  AMREX_ASSERT(!a_pelelm->m_incompressible);
  auto const in_dat = statefab.array();
  auto der = derfab.array(dcomp);
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    der(i, j, k, 0) = in_dat(i, j, k, DENSITY);
    for (int n = 0; n < NUM_SPECIES; n++) {
      der(i, j, k, 0) -= in_dat(i, j, k, FIRSTSPEC + n);
    }
  });
}

//
// Compute cell averaged pressure from nodes
//
void
pelelmex_deravgpress(
  PeleLM* /*a_pelelm*/,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const FArrayBox& /*statefab*/,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& pressfab,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  AMREX_ASSERT(derfab.box().contains(bx));
  auto const in_dat = pressfab.array();
  auto der = derfab.array(dcomp);
  Real factor = 1.0 / (AMREX_D_TERM(2.0, *2.0, *2.0));
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    der(i, j, k) =
      factor * (in_dat(i + 1, j, k) + in_dat(i, j, k)
#if (AMREX_SPACEDIM >= 2)
                + in_dat(i + 1, j + 1, k) + in_dat(i, j + 1, k)
#if (AMREX_SPACEDIM == 3)
                + in_dat(i + 1, j, k + 1) + in_dat(i, j, k + 1) +
                in_dat(i + 1, j + 1, k + 1) + in_dat(i, j + 1, k + 1)
#endif
#endif
               );
  });
}

//
// Compute the velocity magnitude
//
void
pelelmex_dermgvel(
  PeleLM* /*a_pelelm*/,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  auto const vel = statefab.array(VELX);
  auto der = derfab.array(dcomp);
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    der(i, j, k) = std::sqrt((AMREX_D_TERM(
      vel(i, j, k, 0) * vel(i, j, k, 0), +vel(i, j, k, 1) * vel(i, j, k, 1),
      +vel(i, j, k, 2) * vel(i, j, k, 2))));
  });
}

//
// Compute vorticity magnitude
//
void
pelelmex_dermgvort(
  PeleLM* /*a_pelelm*/,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& geom,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);
               , const amrex::Real idy = geom.InvCellSize(1);
               , const amrex::Real idz = geom.InvCellSize(2););

  auto const& dat_arr = statefab.const_array();
  auto const& vort_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
  const auto& ebfab = static_cast<EBFArrayBox const&>(statefab);
  const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

  auto typ = flags.getType(bx);

  if (typ == FabType::covered) {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      vort_arr(i, j, k) = 0.0;
    });
  } else if (typ == FabType::singlevalued) {
    const auto& flag_fab = flags.const_array();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      constexpr amrex::Real c0 = -1.5;
      constexpr amrex::Real c1 = 2.0;
      constexpr amrex::Real c2 = -0.5;
      if (flag_fab(i, j, k).isCovered()) {
        vort_arr(i, j, k) = 0.0;
      } else {
        // Define interpolation lambda
        auto onesided =
          [](const Real& v0, const Real& v1, const Real& v2) -> Real {
          return c0 * v0 + c1 * v1 + c2 * v2;
        };

        amrex::Real vx = 0.0;
        amrex::Real uy = 0.0;
#if (AMREX_SPACEDIM == 2)
        // Need to check if there are covered cells in neighbours --
        // -- if so, use one-sided difference computation (but still quadratic)
        if (!flag_fab(i, j, k).isConnected(1, 0, 0)) {
          vx = -onesided(
                 dat_arr(i, j, k, 1), dat_arr(i - 1, j, k, 1),
                 dat_arr(i - 2, j, k, 1)) *
               idx;
        } else if (!flag_fab(i, j, k).isConnected(-1, 0, 0)) {
          vx = onesided(
                 dat_arr(i, j, k, 1), dat_arr(i + 1, j, k, 1),
                 dat_arr(i + 2, j, k, 1)) *
               idx;
        } else {
          vx = 0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
        }
        // Do the same in y-direction
        if (!flag_fab(i, j, k).isConnected(0, 1, 0)) {
          uy = -onesided(
                 dat_arr(i, j, k, 0), dat_arr(i, j - 1, k, 0),
                 dat_arr(i, j - 2, k, 0)) *
               idy;
        } else if (!flag_fab(i, j, k).isConnected(0, -1, 0)) {
          uy = onesided(
                 dat_arr(i, j, k, 0), dat_arr(i, j + 1, k, 0),
                 dat_arr(i, j + 2, k, 0)) *
               idy;
        } else {
          uy = 0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
        }
        vort_arr(i, j, k) = std::abs(vx - uy);

#elif (AMREX_SPACEDIM == 3)
               amrex::Real wx = 0.0;
               amrex::Real wy = 0.0;
               amrex::Real uz = 0.0;
               amrex::Real vz = 0.0;
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   // Covered cell to the right, go fish left
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
                   wx = - onesided(dat_arr(i,j,k,2), dat_arr(i-1,j,k,2), dat_arr(i-2,j,k,2)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   // Covered cell to the left, go fish right
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
                   wx = onesided(dat_arr(i,j,k,2), dat_arr(i+1,j,k,2), dat_arr(i+2,j,k,2)) * idx;
               } else {
                   // No covered cells right or left, use standard stencil
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                   wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
                   wy = - onesided(dat_arr(i,j,k,2), dat_arr(i,j-1,k,2), dat_arr(i,j-2,k,2)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
                   wy = onesided(dat_arr(i,j,k,2), dat_arr(i,j+1,k,2), dat_arr(i,j+2,k,2)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                   wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
               }
               // Do the same in z-direction
               if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                   uz = - onesided(dat_arr(i,j,k,0), dat_arr(i,j,k-1,0), dat_arr(i,j,k-2,0)) * idz;
                   vz = - onesided(dat_arr(i,j,k,1), dat_arr(i,j,k-1,1), dat_arr(i,j,k-2,1)) * idz;
               } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                   uz = onesided(dat_arr(i,j,k,0), dat_arr(i,j,k+1,0), dat_arr(i,j,k+2,0)) * idz;
                   vz = onesided(dat_arr(i,j,k,1), dat_arr(i,j,k+1,1), dat_arr(i,j,k+2,1)) * idz;
               } else {
                   uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                   vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
               }
               vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
      }
    });
  } else
#endif // Check on EB
  {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if (AMREX_SPACEDIM == 2)
      amrex::Real vx =
        0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
      amrex::Real uy =
        0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
      vort_arr(i, j, k) = std::abs(vx - uy);

#elif (AMREX_SPACEDIM == 3)
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            vort_arr(i,j,k) = std::sqrt((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
    });
  }
}

//
// Compute vorticity components
//
void
pelelmex_dervort(
  PeleLM* /*a_pelelm*/,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& geom,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  amrex::ignore_unused(ncomp, bx);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);
               , const amrex::Real idy = geom.InvCellSize(1);
               , const amrex::Real idz = geom.InvCellSize(2););

  auto const& dat_arr = statefab.const_array();
  auto const& vort_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
  const auto& ebfab = static_cast<EBFArrayBox const&>(statefab);
  const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

  auto typ = flags.getType(bx);

  if (typ == FabType::covered) {
    amrex::ParallelFor(
      bx, ncomp, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
        vort_arr(i, j, k, n) = 0.0;
      });
  } else if (typ == FabType::singlevalued) {
    const auto& flag_fab = flags.const_array();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      constexpr amrex::Real c0 = -1.5;
      constexpr amrex::Real c1 = 2.0;
      constexpr amrex::Real c2 = -0.5;
      if (flag_fab(i, j, k).isCovered()) {
        for (int n{0}; n < ncomp; ++n) {
          vort_arr(i, j, k, n) = 0.0;
        }
      } else {
        // Define interpolation lambda
        auto onesided =
          [](const Real& v0, const Real& v1, const Real& v2) -> Real {
          return c0 * v0 + c1 * v1 + c2 * v2;
        };

        amrex::Real vx = 0.0;
        amrex::Real uy = 0.0;
#if (AMREX_SPACEDIM == 2)
        // Need to check if there are covered cells in neighbours --
        // -- if so, use one-sided difference computation (but still quadratic)
        if (!flag_fab(i, j, k).isConnected(1, 0, 0)) {
          vx = -onesided(
                 dat_arr(i, j, k, 1), dat_arr(i - 1, j, k, 1),
                 dat_arr(i - 2, j, k, 1)) *
               idx;
        } else if (!flag_fab(i, j, k).isConnected(-1, 0, 0)) {
          vx = onesided(
                 dat_arr(i, j, k, 1), dat_arr(i + 1, j, k, 1),
                 dat_arr(i + 2, j, k, 1)) *
               idx;
        } else {
          vx = 0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
        }
        // Do the same in y-direction
        if (!flag_fab(i, j, k).isConnected(0, 1, 0)) {
          uy = -onesided(
                 dat_arr(i, j, k, 0), dat_arr(i, j - 1, k, 0),
                 dat_arr(i, j - 2, k, 0)) *
               idy;
        } else if (!flag_fab(i, j, k).isConnected(0, -1, 0)) {
          uy = onesided(
                 dat_arr(i, j, k, 0), dat_arr(i, j + 1, k, 0),
                 dat_arr(i, j + 2, k, 0)) *
               idy;
        } else {
          uy = 0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
        }
        vort_arr(i, j, k) = vx - uy;

#elif (AMREX_SPACEDIM == 3)
               amrex::Real wx = 0.0;
               amrex::Real wy = 0.0;
               amrex::Real uz = 0.0;
               amrex::Real vz = 0.0;
               // Need to check if there are covered cells in neighbours --
               // -- if so, use one-sided difference computation (but still quadratic)
               if (!flag_fab(i,j,k).isConnected( 1,0,0)) {
                   // Covered cell to the right, go fish left
                   vx = - onesided(dat_arr(i,j,k,1), dat_arr(i-1,j,k,1), dat_arr(i-2,j,k,1)) * idx;
                   wx = - onesided(dat_arr(i,j,k,2), dat_arr(i-1,j,k,2), dat_arr(i-2,j,k,2)) * idx;
               } else if (!flag_fab(i,j,k).isConnected(-1,0,0)) {
                   // Covered cell to the left, go fish right
                   vx = onesided(dat_arr(i,j,k,1), dat_arr(i+1,j,k,1), dat_arr(i+2,j,k,1)) * idx;
                   wx = onesided(dat_arr(i,j,k,2), dat_arr(i+1,j,k,2), dat_arr(i+2,j,k,2)) * idx;
               } else {
                   // No covered cells right or left, use standard stencil
                   vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
                   wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;
               }
               // Do the same in y-direction
               if (!flag_fab(i,j,k).isConnected(0, 1,0)) {
                   uy = - onesided(dat_arr(i,j,k,0), dat_arr(i,j-1,k,0), dat_arr(i,j-2,k,0)) * idy;
                   wy = - onesided(dat_arr(i,j,k,2), dat_arr(i,j-1,k,2), dat_arr(i,j-2,k,2)) * idy;
               } else if (!flag_fab(i,j,k).isConnected(0,-1,0)) {
                   uy = onesided(dat_arr(i,j,k,0), dat_arr(i,j+1,k,0), dat_arr(i,j+2,k,0)) * idy;
                   wy = onesided(dat_arr(i,j,k,2), dat_arr(i,j+1,k,2), dat_arr(i,j+2,k,2)) * idy;
               } else {
                   uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
                   wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;
               }
               // Do the same in z-direction
               if (!flag_fab(i,j,k).isConnected(0,0, 1)) {
                   uz = - onesided(dat_arr(i,j,k,0), dat_arr(i,j,k-1,0), dat_arr(i,j,k-2,0)) * idz;
                   vz = - onesided(dat_arr(i,j,k,1), dat_arr(i,j,k-1,1), dat_arr(i,j,k-2,1)) * idz;
               } else if (!flag_fab(i,j,k).isConnected(0,0,-1)) {
                   uz = onesided(dat_arr(i,j,k,0), dat_arr(i,j,k+1,0), dat_arr(i,j,k+2,0)) * idz;
                   vz = onesided(dat_arr(i,j,k,1), dat_arr(i,j,k+1,1), dat_arr(i,j,k+2,1)) * idz;
               } else {
                   uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
                   vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;
               }
               vort_arr(i,j,k,0) = (wy-vz)*(wy-vz);
               vort_arr(i,j,k,1) = (uz-wx)*(uz-wx);
               vort_arr(i,j,k,2) = (vx-uy)*(vx-uy);
#endif
      }
    });
  } else
#endif // Check on EB

  {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
#if (AMREX_SPACEDIM == 2)
      amrex::Real vx =
        0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
      amrex::Real uy =
        0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
      vort_arr(i, j, k) = vx - uy;

#elif (AMREX_SPACEDIM == 3)
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            vort_arr(i,j,k,0) = (wy-vz)*(wy-vz);
            vort_arr(i,j,k,1) = (uz-wx)*(uz-wx);
            vort_arr(i,j,k,2) = (vx-uy)*(vx-uy);
#endif
    });
  }
}

//
// Compute cell-centered coordinates
//
void
pelelmex_dercoord(
  PeleLM* /*a_pelelm*/,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox&
#ifdef AMREX_USE_EB
    statefab
#else
/*unused*/
#endif
  ,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& geom,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::ignore_unused(ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  AMREX_D_TERM(const amrex::Real dx = geom.CellSize(0);
               , const amrex::Real dy = geom.CellSize(1);
               , const amrex::Real dz = geom.CellSize(2););

  auto const& coord_arr = derfab.array(dcomp);
  const auto geomdata = geom.data();

#ifdef AMREX_USE_EB
  AMREX_ASSERT(statefab.box().contains(bx));
  const auto& ebfab = static_cast<EBFArrayBox const&>(statefab);
  const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

  auto typ = flags.getType(bx);
  // Compute cell center coordinates even in covered boxes/cell. Only
  // modify the cell-center in cut cells
  if (typ == FabType::singlevalued) {
    const auto& flag_arr = flags.const_array();
    const auto& ccent_fab = ebfab.getCentroidData();
    const auto& ccent_arr = ccent_fab->const_array();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      if (flag_arr(i, j, k).isCovered() || flag_arr(i, j, k).isRegular()) {
        AMREX_D_TERM(coord_arr(i, j, k, 0) = prob_lo[0] + (i + 0.5) * dx;
                     , coord_arr(i, j, k, 1) = prob_lo[1] + (j + 0.5) * dy;
                     , coord_arr(i, j, k, 2) = prob_lo[2] + (k + 0.5) * dz;);
      } else {
        AMREX_D_TERM(coord_arr(i, j, k, 0) =
                       prob_lo[0] + (i + 0.5 + ccent_arr(i, j, k, 0)) * dx;
                     , coord_arr(i, j, k, 1) =
                         prob_lo[1] + (j + 0.5 + ccent_arr(i, j, k, 1)) * dy;
                     , coord_arr(i, j, k, 2) =
                         prob_lo[2] + (k + 0.5 + ccent_arr(i, j, k, 2)) * dz;);
      }
    });
  } else
#endif
  {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const amrex::Real* prob_lo = geomdata.ProbLo();
      AMREX_D_TERM(coord_arr(i, j, k, 0) = prob_lo[0] + (i + 0.5) * dx;
                   , coord_arr(i, j, k, 1) = prob_lo[1] + (j + 0.5) * dy;
                   , coord_arr(i, j, k, 2) = prob_lo[2] + (k + 0.5) * dz;);
    });
  }
}

//
// Compute Q-criterion
//
void
pelelmex_derQcrit(
  PeleLM* /*a_pelelm*/,
  const Box&
#if AMREX_SPACEDIM == 3
    bx
#endif
  ,
  FArrayBox&
#if AMREX_SPACEDIM == 3
    derfab
#endif
  ,
  int
#if AMREX_SPACEDIM == 3
    dcomp
#endif
  ,
  int /*ncomp*/,
  const FArrayBox&
#if AMREX_SPACEDIM == 3
    statefab
#endif
  ,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry&
#if AMREX_SPACEDIM == 3
    geom
#endif
  ,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
#if AMREX_SPACEDIM == 3
  AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);
               , const amrex::Real idy = geom.InvCellSize(1);
               , const amrex::Real idz = geom.InvCellSize(2););

  auto const& dat_arr = statefab.const_array();
  auto const& qcrit_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
  const auto& ebfab = static_cast<EBFArrayBox const&>(statefab);
  const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

  auto typ = flags.getType(bx);

  if (typ == FabType::covered) {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      qcrit_arr(i, j, k) = 0.0;
    });
  } else if (typ == FabType::singlevalued) {
    const auto& flag_fab = flags.const_array();
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (flag_fab(i, j, k).isCovered()) {
        qcrit_arr(i, j, k) = 0.0;
      } else {
        // Define interpolation lambda
        constexpr amrex::Real c0 = -1.5;
        constexpr amrex::Real c1 = 2.0;
        constexpr amrex::Real c2 = -0.5;
        auto onesided =
          [](const Real& v0, const Real& v1, const Real& v2) -> Real {
          return c0 * v0 + c1 * v1 + c2 * v2;
        };

        // Strain rate tensor
        Array2D<Real, 0, 2, 0, 2> gradU;
        if (!flag_fab(i, j, k).isConnected(1, 0, 0)) {
          gradU(0, 0) = -onesided(
                          dat_arr(i, j, k, 0), dat_arr(i - 1, j, k, 0),
                          dat_arr(i - 2, j, k, 0)) *
                        idx;
          gradU(1, 0) = -onesided(
                          dat_arr(i, j, k, 1), dat_arr(i - 1, j, k, 1),
                          dat_arr(i - 2, j, k, 1)) *
                        idx;
          gradU(2, 0) = -onesided(
                          dat_arr(i, j, k, 2), dat_arr(i - 1, j, k, 2),
                          dat_arr(i - 2, j, k, 2)) *
                        idx;
        } else if (!flag_fab(i, j, k).isConnected(-1, 0, 0)) {
          gradU(0, 0) = onesided(
                          dat_arr(i, j, k, 0), dat_arr(i + 1, j, k, 0),
                          dat_arr(i + 2, j, k, 0)) *
                        idx;
          gradU(1, 0) = onesided(
                          dat_arr(i, j, k, 1), dat_arr(i + 1, j, k, 1),
                          dat_arr(i + 2, j, k, 1)) *
                        idx;
          gradU(2, 0) = onesided(
                          dat_arr(i, j, k, 2), dat_arr(i + 1, j, k, 2),
                          dat_arr(i + 2, j, k, 2)) *
                        idx;
        } else {
          gradU(0, 0) =
            0.5 * (dat_arr(i + 1, j, k, 0) - dat_arr(i - 1, j, k, 0)) * idx;
          gradU(1, 0) =
            0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
          gradU(2, 0) =
            0.5 * (dat_arr(i + 1, j, k, 2) - dat_arr(i - 1, j, k, 2)) * idx;
        }
        if (!flag_fab(i, j, k).isConnected(0, 1, 0)) {
          gradU(0, 1) = -onesided(
                          dat_arr(i, j, k, 0), dat_arr(i, j - 1, k, 0),
                          dat_arr(i, j - 2, k, 0)) *
                        idy;
          gradU(1, 1) = -onesided(
                          dat_arr(i, j, k, 1), dat_arr(i, j - 1, k, 1),
                          dat_arr(i, j - 2, k, 1)) *
                        idy;
          gradU(2, 1) = -onesided(
                          dat_arr(i, j, k, 2), dat_arr(i, j - 1, k, 2),
                          dat_arr(i, j - 2, k, 2)) *
                        idy;
        } else if (!flag_fab(i, j, k).isConnected(0, -1, 0)) {
          gradU(0, 1) = onesided(
                          dat_arr(i, j, k, 0), dat_arr(i, j + 1, k, 0),
                          dat_arr(i, j + 2, k, 0)) *
                        idy;
          gradU(1, 1) = onesided(
                          dat_arr(i, j, k, 1), dat_arr(i, j + 1, k, 1),
                          dat_arr(i, j + 2, k, 1)) *
                        idy;
          gradU(2, 1) = onesided(
                          dat_arr(i, j, k, 2), dat_arr(i, j + 1, k, 2),
                          dat_arr(i, j + 2, k, 2)) *
                        idy;
        } else {
          gradU(0, 1) =
            0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
          gradU(1, 1) =
            0.5 * (dat_arr(i, j + 1, k, 1) - dat_arr(i, j - 1, k, 1)) * idy;
          gradU(2, 1) =
            0.5 * (dat_arr(i, j + 1, k, 2) - dat_arr(i, j - 1, k, 2)) * idy;
        }
        if (!flag_fab(i, j, k).isConnected(0, 0, 1)) {
          gradU(0, 2) = -onesided(
                          dat_arr(i, j, k, 0), dat_arr(i, j, k - 1, 0),
                          dat_arr(i, j, k - 2, 0)) *
                        idz;
          gradU(1, 2) = -onesided(
                          dat_arr(i, j, k, 1), dat_arr(i, j, k - 1, 1),
                          dat_arr(i, j, k - 2, 1)) *
                        idz;
          gradU(2, 2) = -onesided(
                          dat_arr(i, j, k, 2), dat_arr(i, j, k - 1, 2),
                          dat_arr(i, j, k - 2, 2)) *
                        idz;
        } else if (!flag_fab(i, j, k).isConnected(0, 0, -1)) {
          gradU(0, 2) = onesided(
                          dat_arr(i, j, k, 0), dat_arr(i, j, k + 1, 0),
                          dat_arr(i, j, k + 2, 0)) *
                        idz;
          gradU(1, 2) = onesided(
                          dat_arr(i, j, k, 1), dat_arr(i, j, k + 1, 1),
                          dat_arr(i, j, k + 2, 1)) *
                        idz;
          gradU(2, 2) = onesided(
                          dat_arr(i, j, k, 2), dat_arr(i, j, k + 1, 2),
                          dat_arr(i, j, k + 2, 2)) *
                        idz;
        } else {
          gradU(0, 2) =
            0.5 * (dat_arr(i, j, k + 1, 0) - dat_arr(i, j, k - 1, 0)) * idz;
          gradU(1, 2) =
            0.5 * (dat_arr(i, j, k + 1, 1) - dat_arr(i, j, k - 1, 1)) * idz;
          gradU(2, 2) =
            0.5 * (dat_arr(i, j, k + 1, 2) - dat_arr(i, j, k - 1, 2)) * idz;
        }

        // Divu
        amrex::Real divU = gradU(0, 0) + gradU(1, 1) + gradU(2, 2);

        // Directly Assemble Sym. & AntiSym. into Qcrit.
        // Remove divU (dilatation) from the Sym. tensor (due to mixing/reaction
        // most often)
        qcrit_arr(i, j, k) = 0.0;
        for (int dim1 = 0; dim1 < AMREX_SPACEDIM; ++dim1) {
          for (int dim2 = 0; dim2 < AMREX_SPACEDIM; ++dim2) {
            Real Ohm = 0.5 * (gradU(dim1, dim2) - gradU(dim2, dim1));
            Real Sij = 0.5 * (gradU(dim1, dim2) + gradU(dim2, dim1));
            if (dim1 == dim2) {
              Sij -= divU / AMREX_SPACEDIM;
            }
            qcrit_arr(i, j, k) += Ohm * Ohm - Sij * Sij;
          }
        }
      }
    });
  } else
#endif
  {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      // Strain rate tensor
      Array2D<Real, 0, 2, 0, 2> gradU;
      gradU(0, 0) =
        0.5 * (dat_arr(i + 1, j, k, 0) - dat_arr(i - 1, j, k, 0)) * idx;
      gradU(0, 1) =
        0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
      gradU(0, 2) =
        0.5 * (dat_arr(i, j, k + 1, 0) - dat_arr(i, j, k - 1, 0)) * idz;
      gradU(1, 0) =
        0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
      gradU(1, 1) =
        0.5 * (dat_arr(i, j + 1, k, 1) - dat_arr(i, j - 1, k, 1)) * idy;
      gradU(1, 2) =
        0.5 * (dat_arr(i, j, k + 1, 1) - dat_arr(i, j, k - 1, 1)) * idz;
      gradU(2, 0) =
        0.5 * (dat_arr(i + 1, j, k, 2) - dat_arr(i - 1, j, k, 2)) * idx;
      gradU(2, 1) =
        0.5 * (dat_arr(i, j + 1, k, 2) - dat_arr(i, j - 1, k, 2)) * idy;
      gradU(2, 2) =
        0.5 * (dat_arr(i, j, k + 1, 2) - dat_arr(i, j, k - 1, 2)) * idz;

      // Divu
      amrex::Real divU = gradU(0, 0) + gradU(1, 1) + gradU(2, 2);

      // Directly Assemble Sym. & AntiSym. into Qcrit.
      // Remove divU (dilatation) from the Sym. tensor (due to mixing/reaction
      // most often)
      qcrit_arr(i, j, k) = 0.0;
      for (int dim1 = 0; dim1 < AMREX_SPACEDIM; ++dim1) {
        for (int dim2 = 0; dim2 < AMREX_SPACEDIM; ++dim2) {
          Real Ohm = 0.5 * (gradU(dim1, dim2) - gradU(dim2, dim1));
          Real Sij = 0.5 * (gradU(dim1, dim2) + gradU(dim2, dim1));
          if (dim1 == dim2) {
            Sij -= divU / AMREX_SPACEDIM;
          }
          qcrit_arr(i, j, k) += Ohm * Ohm - Sij * Sij;
        }
      }
    });
  }
#endif
}

//
// Compute the kinetic energy
//
void
pelelmex_derkineticenergy(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  if (a_pelelm->m_incompressible != 0) {
    auto const vel = statefab.array(VELX);
    auto der = derfab.array(dcomp);
    amrex::ParallelFor(
      bx, [=, rho = a_pelelm->m_rho] AMREX_GPU_DEVICE(
            int i, int j, int k) noexcept {
        der(i, j, k) = 0.5 * rho *
                       (AMREX_D_TERM(
                         vel(i, j, k, 0) * vel(i, j, k, 0),
                         +vel(i, j, k, 1) * vel(i, j, k, 1),
                         +vel(i, j, k, 2) * vel(i, j, k, 2)));
      });
  } else {
    auto const rho = statefab.array(DENSITY);
    auto const vel = statefab.array(VELX);
    auto der = derfab.array(dcomp);
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      der(i, j, k) =
        0.5 * rho(i, j, k) *
        (AMREX_D_TERM(
          vel(i, j, k, 0) * vel(i, j, k, 0), +vel(i, j, k, 1) * vel(i, j, k, 1),
          +vel(i, j, k, 2) * vel(i, j, k, 2)));
    });
  }
}

//
// Compute enstrophy
//
void
pelelmex_derenstrophy(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& geom,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  AMREX_D_TERM(const amrex::Real idx = geom.InvCellSize(0);
               , const amrex::Real idy = geom.InvCellSize(1);
               , const amrex::Real idz = geom.InvCellSize(2););

  auto const& dat_arr = statefab.const_array(VELX);
  auto const& rho_arr = (a_pelelm->m_incompressible) != 0
                          ? Array4<const Real>{}
                          : statefab.const_array(DENSITY);
  auto const& ens_arr = derfab.array(dcomp);

#ifdef AMREX_USE_EB
  const auto& ebfab = static_cast<EBFArrayBox const&>(statefab);
  const EBCellFlagFab& flags = ebfab.getEBCellFlagFab();

  auto typ = flags.getType(bx);

  if (typ == FabType::covered) {
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      ens_arr(i, j, k) = 0.0;
    });
  } else if (typ == FabType::singlevalued) {
    const auto& flag_fab = flags.const_array();
    amrex::ParallelFor(
      bx,
      [=, incomp = a_pelelm->m_incompressible,
       rho = a_pelelm->m_rho] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        constexpr amrex::Real c0 = -1.5;
        constexpr amrex::Real c1 = 2.0;
        constexpr amrex::Real c2 = -0.5;
        if (flag_fab(i, j, k).isCovered()) {
          ens_arr(i, j, k) = 0.0;
        } else {
          Real l_rho = rho;
          if (incomp == 0) {
            l_rho = rho_arr(i, j, k);
          }
          // Define interpolation lambda
          auto onesided =
            [](const Real& v0, const Real& v1, const Real& v2) -> Real {
            return c0 * v0 + c1 * v1 + c2 * v2;
          };

          amrex::Real vx = 0.0;
          amrex::Real uy = 0.0;
#if (AMREX_SPACEDIM == 2)
          // Need to check if there are covered cells in neighbours --
          // -- if so, use one-sided difference computation (but still
          // quadratic)
          if (!flag_fab(i, j, k).isConnected(1, 0, 0)) {
            vx = -onesided(
                   dat_arr(i, j, k, 1), dat_arr(i - 1, j, k, 1),
                   dat_arr(i - 2, j, k, 1)) *
                 idx;
          } else if (!flag_fab(i, j, k).isConnected(-1, 0, 0)) {
            vx = onesided(
                   dat_arr(i, j, k, 1), dat_arr(i + 1, j, k, 1),
                   dat_arr(i + 2, j, k, 1)) *
                 idx;
          } else {
            vx =
              0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
          }
          // Do the same in y-direction
          if (!flag_fab(i, j, k).isConnected(0, 1, 0)) {
            uy = -onesided(
                   dat_arr(i, j, k, 0), dat_arr(i, j - 1, k, 0),
                   dat_arr(i, j - 2, k, 0)) *
                 idy;
          } else if (!flag_fab(i, j, k).isConnected(0, -1, 0)) {
            uy = onesided(
                   dat_arr(i, j, k, 0), dat_arr(i, j + 1, k, 0),
                   dat_arr(i, j + 2, k, 0)) *
                 idy;
          } else {
            uy =
              0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
          }
          ens_arr(i, j, k) = 0.5 * l_rho * (vx - uy) * (vx - uy);

#elif (AMREX_SPACEDIM == 3)
          amrex::Real wx = 0.0;
          amrex::Real wy = 0.0;
          amrex::Real uz = 0.0;
          amrex::Real vz = 0.0;
          // Need to check if there are covered cells in neighbours --
          // -- if so, use one-sided difference computation (but still
          // quadratic)
          if (!flag_fab(i, j, k).isConnected(1, 0, 0)) {
            // Covered cell to the right, go fish left
            vx = -onesided(
                   dat_arr(i, j, k, 1), dat_arr(i - 1, j, k, 1),
                   dat_arr(i - 2, j, k, 1)) *
                 idx;
            wx = -onesided(
                   dat_arr(i, j, k, 2), dat_arr(i - 1, j, k, 2),
                   dat_arr(i - 2, j, k, 2)) *
                 idx;
          } else if (!flag_fab(i, j, k).isConnected(-1, 0, 0)) {
            // Covered cell to the left, go fish right
            vx = onesided(
                   dat_arr(i, j, k, 1), dat_arr(i + 1, j, k, 1),
                   dat_arr(i + 2, j, k, 1)) *
                 idx;
            wx = onesided(
                   dat_arr(i, j, k, 2), dat_arr(i + 1, j, k, 2),
                   dat_arr(i + 2, j, k, 2)) *
                 idx;
          } else {
            // No covered cells right or left, use standard stencil
            vx =
              0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
            wx =
              0.5 * (dat_arr(i + 1, j, k, 2) - dat_arr(i - 1, j, k, 2)) * idx;
          }
          // Do the same in y-direction
          if (!flag_fab(i, j, k).isConnected(0, 1, 0)) {
            uy = -onesided(
                   dat_arr(i, j, k, 0), dat_arr(i, j - 1, k, 0),
                   dat_arr(i, j - 2, k, 0)) *
                 idy;
            wy = -onesided(
                   dat_arr(i, j, k, 2), dat_arr(i, j - 1, k, 2),
                   dat_arr(i, j - 2, k, 2)) *
                 idy;
          } else if (!flag_fab(i, j, k).isConnected(0, -1, 0)) {
            uy = onesided(
                   dat_arr(i, j, k, 0), dat_arr(i, j + 1, k, 0),
                   dat_arr(i, j + 2, k, 0)) *
                 idy;
            wy = onesided(
                   dat_arr(i, j, k, 2), dat_arr(i, j + 1, k, 2),
                   dat_arr(i, j + 2, k, 2)) *
                 idy;
          } else {
            uy =
              0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
            wy =
              0.5 * (dat_arr(i, j + 1, k, 2) - dat_arr(i, j - 1, k, 2)) * idy;
          }
          // Do the same in z-direction
          if (!flag_fab(i, j, k).isConnected(0, 0, 1)) {
            uz = -onesided(
                   dat_arr(i, j, k, 0), dat_arr(i, j, k - 1, 0),
                   dat_arr(i, j, k - 2, 0)) *
                 idz;
            vz = -onesided(
                   dat_arr(i, j, k, 1), dat_arr(i, j, k - 1, 1),
                   dat_arr(i, j, k - 2, 1)) *
                 idz;
          } else if (!flag_fab(i, j, k).isConnected(0, 0, -1)) {
            uz = onesided(
                   dat_arr(i, j, k, 0), dat_arr(i, j, k + 1, 0),
                   dat_arr(i, j, k + 2, 0)) *
                 idz;
            vz = onesided(
                   dat_arr(i, j, k, 1), dat_arr(i, j, k + 1, 1),
                   dat_arr(i, j, k + 2, 1)) *
                 idz;
          } else {
            uz =
              0.5 * (dat_arr(i, j, k + 1, 0) - dat_arr(i, j, k - 1, 0)) * idz;
            vz =
              0.5 * (dat_arr(i, j, k + 1, 1) - dat_arr(i, j, k - 1, 1)) * idz;
          }
          ens_arr(i, j, k) = 0.5 * l_rho *
                             ((wy - vz) * (wy - vz) + (uz - wx) * (uz - wx) +
                              (vx - uy) * (vx - uy));
#endif
        }
      });
  } else
#endif
  {
    amrex::ParallelFor(
      bx,
      [=, incomp = a_pelelm->m_incompressible,
       rho = a_pelelm->m_rho] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        Real l_rho = rho;
        if (incomp == 0) {
          l_rho = rho_arr(i, j, k);
        }
#if (AMREX_SPACEDIM == 2)
        amrex::Real vx =
          0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
        amrex::Real uy =
          0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
        ens_arr(i, j, k) = 0.5 * l_rho * (vx - uy) * (vx - uy);

#elif (AMREX_SPACEDIM == 3)
        amrex::Real vx =
          0.5 * (dat_arr(i + 1, j, k, 1) - dat_arr(i - 1, j, k, 1)) * idx;
        amrex::Real wx =
          0.5 * (dat_arr(i + 1, j, k, 2) - dat_arr(i - 1, j, k, 2)) * idx;

        amrex::Real uy =
          0.5 * (dat_arr(i, j + 1, k, 0) - dat_arr(i, j - 1, k, 0)) * idy;
        amrex::Real wy =
          0.5 * (dat_arr(i, j + 1, k, 2) - dat_arr(i, j - 1, k, 2)) * idy;

        amrex::Real uz =
          0.5 * (dat_arr(i, j, k + 1, 0) - dat_arr(i, j, k - 1, 0)) * idz;
        amrex::Real vz =
          0.5 * (dat_arr(i, j, k + 1, 1) - dat_arr(i, j, k - 1, 1)) * idz;

        ens_arr(i, j, k) = 0.5 * l_rho *
                           ((wy - vz) * (wy - vz) + (uz - wx) * (uz - wx) +
                            (vx - uy) * (vx - uy));
#endif
      });
  }
}

//
// Compute mixture fraction
//
void
pelelmex_dermixfrac(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::ignore_unused(ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(ncomp == 1);

  if (a_pelelm->Zfu < 0.0) {
    amrex::Abort("Mixture fraction not initialized");
  }

  auto const density = statefab.array(DENSITY);
  auto const rhoY = statefab.array(FIRSTSPEC);
  auto mixt_frac = derfab.array(dcomp);

  amrex::Real Zox_lcl = a_pelelm->Zox;
  amrex::Real Zfu_lcl = a_pelelm->Zfu;
  amrex::Real denom_inv = 1.0 / (Zfu_lcl - Zox_lcl);
  amrex::GpuArray<amrex::Real, NUM_SPECIES> fact_Bilger;
  for (int n = 0; n < NUM_SPECIES; ++n) {
    fact_Bilger[n] = a_pelelm->spec_Bilger_fact[n];
  }

  amrex::ParallelFor(
    bx, [density, rhoY, mixt_frac, fact_Bilger, Zox_lcl,
         denom_inv] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real rho_inv = 1.0_rt / density(i, j, k);
      mixt_frac(i, j, k) = 0.0_rt;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        mixt_frac(i, j, k) += (rhoY(i, j, k, n) * fact_Bilger[n]) * rho_inv;
      }
      mixt_frac(i, j, k) = (mixt_frac(i, j, k) - Zox_lcl) * denom_inv;
    });
}

//
// Compute progress variable
//
void
pelelmex_derprogvar(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  amrex::ignore_unused(ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(ncomp == 1);

  if (a_pelelm->m_C0 < 0.0) {
    amrex::Abort("Progress variable not initialized");
  }

  auto const density = statefab.array(DENSITY);
  auto const rhoY = statefab.array(FIRSTSPEC);
  auto const temp = statefab.array(TEMP);
  auto prog_var = derfab.array(dcomp);

  amrex::Real C0_lcl = a_pelelm->m_C0;
  amrex::Real C1_lcl = a_pelelm->m_C1;
  amrex::Real denom_inv = 1.0 / (C1_lcl - C0_lcl);
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 1> Cweights;
  for (int n = 0; n < NUM_SPECIES + 1; ++n) {
    Cweights[n] = a_pelelm->m_Cweights[n];
  }

  amrex::ParallelFor(
    bx, [=, revert = a_pelelm->m_Crevert] AMREX_GPU_DEVICE(
          int i, int j, int k) noexcept {
      amrex::Real rho_inv = 1.0_rt / density(i, j, k);
      prog_var(i, j, k) = 0.0_rt;
      for (int n = 0; n < NUM_SPECIES; ++n) {
        prog_var(i, j, k) += (rhoY(i, j, k, n) * Cweights[n]) * rho_inv;
      }
      prog_var(i, j, k) += temp(i, j, k) * Cweights[NUM_SPECIES];
      if (revert != 0) {
        prog_var(i, j, k) = 1.0 - (prog_var(i, j, k) - C0_lcl) * denom_inv;
      } else {
        prog_var(i, j, k) = (prog_var(i, j, k) - C0_lcl) * denom_inv;
      }
    });
}

//
// Extract mixture viscosity
//
void
pelelmex_dervisc(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::ignore_unused(ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);

  if (a_pelelm->m_incompressible != 0) {
    derfab.setVal<RunOn::Device>(a_pelelm->m_mu, bx, dcomp, 1);
  } else {
    auto const& rhoY = statefab.const_array(FIRSTSPEC);
    auto const& T = statefab.array(TEMP);
    auto der = derfab.array(dcomp);
    auto const* ltransparm = a_pelelm->trans_parms.device_parm();
    amrex::ParallelFor(
      bx, [rhoY, T, der,
           ltransparm] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        getVelViscosity(i, j, k, rhoY, T, der, ltransparm);
      });
  }
}

//
// Extract mixture averaged species diffusion coefficients
//
void
pelelmex_derdiffc(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::ignore_unused(ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  if (a_pelelm->m_use_soret != 0) {
    AMREX_ASSERT(ncomp == 2 * NUM_SPECIES);
  }
  if (a_pelelm->m_use_soret == 0) {
    AMREX_ASSERT(ncomp == NUM_SPECIES);
  }
  bool do_fixed_Le = (a_pelelm->m_fixed_Le != 0);
  bool do_fixed_Pr = (a_pelelm->m_fixed_Pr != 0);
  bool do_soret = (a_pelelm->m_use_soret != 0);
  FArrayBox dummies(bx, NUM_SPECIES + 2, The_Async_Arena());
  auto const& rhoY = statefab.const_array(FIRSTSPEC);
  auto const& T = statefab.array(TEMP);
  auto rhoD = derfab.array(dcomp);
  auto lambda = dummies.array(0);
  auto mu = dummies.array(1);
  auto const* ltransparm = a_pelelm->trans_parms.device_parm();
  auto const* leosparm = a_pelelm->eos_parms.device_parm();
  auto rhotheta = do_soret ? derfab.array(dcomp + NUM_SPECIES)
                           : dummies.array(2); // dummy for no soret
  amrex::Real LeInv = a_pelelm->m_Lewis_inv;
  amrex::Real PrInv = a_pelelm->m_Prandtl_inv;
  amrex::ParallelFor(
    bx, [do_fixed_Le, do_fixed_Pr, do_soret, LeInv, PrInv, rhoY, T, rhoD,
         rhotheta, lambda, mu, ltransparm,
         leosparm] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      getTransportCoeff(
        i, j, k, do_fixed_Le, do_fixed_Pr, do_soret, LeInv, PrInv, rhoY, T,
        rhoD, rhotheta, lambda, mu, ltransparm, leosparm);
    });
}

//
// Extract thermal diffusivity
//
void
pelelmex_derlambda(
  PeleLM* a_pelelm,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int ncomp,
  const FArrayBox& statefab,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)
{
  amrex::ignore_unused(ncomp);
  AMREX_ASSERT(derfab.box().contains(bx));
  AMREX_ASSERT(statefab.box().contains(bx));
  AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
  bool do_fixed_Le = (a_pelelm->m_fixed_Le != 0);
  bool do_fixed_Pr = (a_pelelm->m_fixed_Pr != 0);
  bool do_soret = (a_pelelm->m_use_soret != 0);
  FArrayBox dummies(bx, 2 * NUM_SPECIES + 1, The_Async_Arena());
  auto const& rhoY = statefab.const_array(FIRSTSPEC);
  auto const& T = statefab.array(TEMP);
  auto rhoD = dummies.array(1);
  auto lambda = derfab.array(dcomp);
  auto mu = dummies.array(0);
  auto rhotheta = dummies.array(NUM_SPECIES + 1);
  auto const* ltransparm = a_pelelm->trans_parms.device_parm();
  auto const* leosparm = a_pelelm->eos_parms.device_parm();
  amrex::Real LeInv = a_pelelm->m_Lewis_inv;
  amrex::Real PrInv = a_pelelm->m_Prandtl_inv;
  amrex::ParallelFor(
    bx, [do_fixed_Le, do_fixed_Pr, do_soret, LeInv, PrInv, rhoY, T, rhoD,
         rhotheta, lambda, mu, ltransparm,
         leosparm] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      getTransportCoeff(
        i, j, k, do_fixed_Le, do_fixed_Pr, do_soret, LeInv, PrInv, rhoY, T,
        rhoD, rhotheta, lambda, mu, ltransparm, leosparm);
    });
}

//
// Extract Distribution Mapping
//
void
pelelmex_derdmap(
  PeleLM* /*a_pelelm*/,
  const Box& bx,
  FArrayBox& derfab,
  int dcomp,
  int /*ncomp*/,
  const FArrayBox& /*statefab*/,
  const FArrayBox& /*reactfab*/,
  const FArrayBox& /*pressfab*/,
  const Geometry& /*geom*/,
  Real /*time*/,
  const Vector<BCRec>& /*bcrec*/,
  int /*level*/)

{
  AMREX_ASSERT(derfab.box().contains(bx));
  auto der = derfab.array(dcomp);
  const int myrank = ParallelDescriptor::MyProc();
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    der(i, j, k) = myrank;
  });
}
