#include <PeleLMDeriveFunc.H>
#include <PeleLM_Index.H>
#include <EOS.H>
#include <mechanism.h>

using namespace amrex;

//
// Extract temp
//
void pelelm_dertemp (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geomdata*/,
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
                         const FArrayBox& statefab, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geomdata*/,
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

//
// Compute cell averaged pressure from nodes
//
void pelelm_deravgpress (const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& /*statefab*/, const FArrayBox& pressfab,
                         const Geometry& /*geomdata*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    auto const in_dat = pressfab.array();
    auto       der = derfab.array(dcomp);
    Real factor = 1.0 / ( AMREX_D_TERM(2.0,*2.0,*2.0) );
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {   
        der(i,j,k) =  factor * (  in_dat(i+1,j,k)     + in_dat(i,j,k)
#if (AMREX_SPACEDIM >= 2 )
                                + in_dat(i+1,j+1,k)   + in_dat(i,j+1,k)
#if (AMREX_SPACEDIM == 3 )
                                + in_dat(i+1,j,k+1)   + in_dat(i,j,k+1)
                                + in_dat(i+1,j+1,k+1) + in_dat(i,j+1,k+1)
#endif
#endif
                                );  
    }); 
}
