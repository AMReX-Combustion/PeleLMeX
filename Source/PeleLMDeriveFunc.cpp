#include <PeleLMDeriveFunc.H>
#include <PeleLM_Index.H>
#include <EOS.H>
#include <mechanism.h>

using namespace amrex;

//
// Extract temp
//
void pelelm_dertemp (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const Geometry& /*geomdata*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        der(i,j,k) = in_dat(i,j,k,TEMP);
    }); 
}

//
// Extract species mass fractions Y_n
//
void pelelm_dermassfrac (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const Geometry& /*geomdata*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx, NUM_SPECIES,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
    {   
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,DENSITY);
        der(i,j,k,n) = in_dat(i,j,k,FIRSTSPEC+n) * rhoinv;
    }); 
}
