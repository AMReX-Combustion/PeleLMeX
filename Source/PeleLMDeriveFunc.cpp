#include <PeleLMDeriveFunc.H>
#include <PeleLM_Index.H>
#include <PelePhysics.H>
#include <mechanism.H>
#include <PeleLM.H>
#include <PeleLM_K.H>

using namespace amrex;

//
// Extract temp
//
void pelelm_dertemp (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geomdata*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        der(i,j,k) = in_dat(i,j,k,TEMP);
    });
}

//
// Compute heat release
//
void pelelm_derheatrelease (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                            const FArrayBox& statefab, const FArrayBox& reactfab, const FArrayBox& /*pressfab*/,
                            const Geometry& /*geomdata*/,
                            Real /*time*/, const Vector<BCRec>& /*bcrec*/, int level)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(!a_pelelm->m_incompressible);

    FArrayBox EnthFab;
    EnthFab.resize(bx,NUM_SPECIES,The_Async_Arena());

    auto const temp = statefab.const_array(TEMP);
    auto const react = reactfab.const_array(0);
    auto const& Hi    = EnthFab.array();
    auto       HRR = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        getHGivenT( i, j, k, temp, Hi );
        HRR(i,j,k) = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
           HRR(i,j,k) -= Hi(i,j,k,n) * react(i,j,k,n);
        }
    });
}

//
// Extract species mass fractions Y_n
//
void pelelm_dermassfrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geomdata*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
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
// Extract species mole fractions X_n
//
void pelelm_dermolefrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                         const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                         const Geometry& /*geomdata*/,
                         Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(statefab.nComp() >= NUM_SPECIES+1);
    AMREX_ASSERT(ncomp == NUM_SPECIES);
    AMREX_ASSERT(!a_pelelm->m_incompressible);
    auto const in_dat = statefab.array();
    auto       der = derfab.array(dcomp);
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real Yt[NUM_SPECIES] = {0.0};
        amrex::Real Xt[NUM_SPECIES] = {0.0};
        amrex::Real rhoinv = 1.0 / in_dat(i,j,k,DENSITY);
        for (int n = 0; n < NUM_SPECIES; n++) {
           Yt[n] = in_dat(i,j,k,FIRSTSPEC+n) * rhoinv;
        }
        auto eos = pele::physics::PhysicsType::eos();
        eos.Y2X(Yt,Xt);
        for (int n = 0; n < NUM_SPECIES; n++) {
           der(i,j,k,n) = Xt[n];
        }
    });
}

//
// Compute cell averaged pressure from nodes
//
void pelelm_deravgpress (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                         const FArrayBox& /*statefab*/, const FArrayBox& /*reactfab*/, const FArrayBox& pressfab,
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

//
// Compute vorticity magnitude
//
void pelelm_dermgvort (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                       const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                       const Geometry& geomdata,
                       Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{

    AMREX_D_TERM(const amrex::Real idx = geomdata.InvCellSize(0);,
                 const amrex::Real idy = geomdata.InvCellSize(1);,
                 const amrex::Real idz = geomdata.InvCellSize(2););

    auto const& dat_arr = statefab.const_array();
    auto const&vort_arr = derfab.array(dcomp);

    // TODO : EB
    // TODO : BCs

    {
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            vort_arr(i,j,k) = vx-uy;

#elif ( AMREX_SPACEDIM == 3 )
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
// Compute the kinetic energy
//
void pelelm_derkineticenergy (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                              const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                              const Geometry& /*geomdata*/,
                              Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    if (a_pelelm->m_incompressible) {
       auto const vel = statefab.array(VELX);
       auto       der = derfab.array(dcomp);
       amrex::ParallelFor(bx,
       [=,rho=a_pelelm->m_rho] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           der(i,j,k) = 0.5 * rho
                            * ( AMREX_D_TERM(  vel(i,j,k,0)*vel(i,j,k,0),
                                             + vel(i,j,k,1)*vel(i,j,k,1),
                                             + vel(i,j,k,2)*vel(i,j,k,2)) );
       });
    } else {
       auto const rho = statefab.array(DENSITY);
       auto const vel = statefab.array(VELX);
       auto       der = derfab.array(dcomp);
       amrex::ParallelFor(bx,
       [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
       {
           der(i,j,k) = 0.5 * rho(i,j,k)
                            * ( AMREX_D_TERM(  vel(i,j,k,0)*vel(i,j,k,0),
                                             + vel(i,j,k,1)*vel(i,j,k,1),
                                             + vel(i,j,k,2)*vel(i,j,k,2)) );
       });
    }
}

//
// Compute enstrophy
//
void pelelm_derenstrophy (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int /*ncomp*/,
                          const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                          const Geometry& geomdata,
                          Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{

    AMREX_D_TERM(const amrex::Real idx = geomdata.InvCellSize(0);,
                 const amrex::Real idy = geomdata.InvCellSize(1);,
                 const amrex::Real idz = geomdata.InvCellSize(2););

    // TODO : EB
    // TODO : BCs

    if (a_pelelm->m_incompressible) {
        auto const&  dat_arr = statefab.const_array(VELX);
        auto const&  ens_arr = derfab.array(dcomp);
        amrex::ParallelFor(bx, [=,rho=a_pelelm->m_rho] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            ens_arr(i,j,k) = 0.5 * rho * (vx-uy)*(vx-uy);

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            ens_arr(i,j,k) = 0.5 * rho * ((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    } else {
        auto const&  dat_arr = statefab.const_array(VELX);
        auto const&  rho_arr = statefab.const_array(DENSITY);
        auto const&  ens_arr = derfab.array(dcomp);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
#if ( AMREX_SPACEDIM == 2 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            ens_arr(i,j,k) = 0.5 * rho_arr(i,j,k) * (vx-uy)*(vx-uy);

#elif ( AMREX_SPACEDIM == 3 )
            amrex::Real vx = 0.5 * (dat_arr(i+1,j,k,1) - dat_arr(i-1,j,k,1)) * idx;
            amrex::Real wx = 0.5 * (dat_arr(i+1,j,k,2) - dat_arr(i-1,j,k,2)) * idx;

            amrex::Real uy = 0.5 * (dat_arr(i,j+1,k,0) - dat_arr(i,j-1,k,0)) * idy;
            amrex::Real wy = 0.5 * (dat_arr(i,j+1,k,2) - dat_arr(i,j-1,k,2)) * idy;

            amrex::Real uz = 0.5 * (dat_arr(i,j,k+1,0) - dat_arr(i,j,k-1,0)) * idz;
            amrex::Real vz = 0.5 * (dat_arr(i,j,k+1,1) - dat_arr(i,j,k-1,1)) * idz;

            ens_arr(i,j,k) = 0.5 * rho_arr(i,j,k) * ((wy-vz)*(wy-vz) + (uz-wx)*(uz-wx) + (vx-uy)*(vx-uy));
#endif
        });
    }

}


//
// Compute mixture fraction
//
void pelelm_dermixfrac (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geomdata*/,
                        Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(ncomp == 1);

    if (a_pelelm->Zfu < 0.0) amrex::Abort("Mixture fraction not initialized");

    auto const density   = statefab.array(DENSITY);
    auto const rhoY      = statefab.array(FIRSTSPEC);
    auto       mixt_frac = derfab.array(dcomp);

    amrex::Real Zox_lcl   = a_pelelm->Zox;
    amrex::Real Zfu_lcl   = a_pelelm->Zfu;
    amrex::Real denom_inv = 1.0 / (Zfu_lcl - Zox_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES> fact_Bilger;
    for (int n=0; n<NUM_SPECIES; ++n) {
        fact_Bilger[n] = a_pelelm->spec_Bilger_fact[n];
    }

    amrex::ParallelFor(bx,
    [density, rhoY, mixt_frac, fact_Bilger, Zox_lcl, denom_inv] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real rho_inv = 1.0_rt / density(i,j,k);
        mixt_frac(i,j,k) = 0.0_rt;
        for (int n = 0; n<NUM_SPECIES; ++n) {
            mixt_frac(i,j,k) += ( rhoY(i,j,k,n) * fact_Bilger[n] ) * rho_inv;
        }
        mixt_frac(i,j,k) = ( mixt_frac(i,j,k) - Zox_lcl ) * denom_inv;
    });
}

//
// Compute progress variable
//
void pelelm_derprogvar (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                        const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                        const Geometry& /*geomdata*/,
                        Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)

{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(ncomp == 1);

    if (a_pelelm->m_C0 < 0.0) amrex::Abort("Progress variable not initialized");

    auto const density  = statefab.array(DENSITY);
    auto const rhoY     = statefab.array(FIRSTSPEC);
    auto const temp     = statefab.array(TEMP);
    auto       prog_var = derfab.array(dcomp);

    amrex::Real C0_lcl   = a_pelelm->m_C0;
    amrex::Real C1_lcl   = a_pelelm->m_C1;
    amrex::Real denom_inv = 1.0 / (C1_lcl - C0_lcl);
    amrex::GpuArray<amrex::Real,NUM_SPECIES+1> Cweights;
    for (int n=0; n<NUM_SPECIES+1; ++n) {
        Cweights[n] = a_pelelm->m_Cweights[n];
    }

    amrex::ParallelFor(bx, [=,revert=a_pelelm->m_Crevert]
    AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        amrex::Real rho_inv = 1.0_rt / density(i,j,k);
        prog_var(i,j,k) = 0.0_rt;
        for (int n = 0; n<NUM_SPECIES; ++n) {
            prog_var(i,j,k) += ( rhoY(i,j,k,n) * Cweights[n] ) * rho_inv;
        }
        prog_var(i,j,k) += temp(i,j,k) * Cweights[NUM_SPECIES];
        if (revert) {
           prog_var(i,j,k) = 1.0 - ( prog_var(i,j,k) - C0_lcl ) * denom_inv;
        } else {
           prog_var(i,j,k) = ( prog_var(i,j,k) - C0_lcl ) * denom_inv;
        }
    });
}

//
// Extract mixture viscosity
//
void pelelm_dervisc (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                     const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                     const Geometry& /*geomdata*/,
                     Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);

    if (a_pelelm->m_incompressible) {
        derfab.setVal<RunOn::Device>(a_pelelm->m_mu,bx,dcomp,1);
    } else {
        auto const& rhoY = statefab.const_array(FIRSTSPEC);
        auto const& T    = statefab.array(TEMP);
        auto       der   = derfab.array(dcomp);
        auto const* ltransparm = a_pelelm->trans_parms.device_trans_parm();
        amrex::ParallelFor(bx,
        [rhoY,T,der,ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            getVelViscosity(i, j, k, rhoY, T, der, ltransparm);
        });
    }
}

//
// Extract mixture averaged species diffusion coefficients
//
void pelelm_derdiffc (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                      const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                      const Geometry& /*geomdata*/,
                      Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);
    AMREX_ASSERT(ncomp == NUM_SPECIES);

    FArrayBox dummies(bx,2,The_Async_Arena());
    auto const& rhoY = statefab.const_array(FIRSTSPEC);
    auto const& T    = statefab.array(TEMP);
    auto       rhoD  = derfab.array(dcomp);
    auto     lambda  = dummies.array(0);
    auto         mu  = dummies.array(1);
    auto const* ltransparm = a_pelelm->trans_parms.device_trans_parm();
    amrex::ParallelFor(bx,
    [rhoY,T,rhoD,lambda,mu,ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        getTransportCoeff(i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
    });
}

//
// Extract thermal diffusivity
//
void pelelm_derlambda (PeleLM* a_pelelm, const Box& bx, FArrayBox& derfab, int dcomp, int ncomp,
                       const FArrayBox& statefab, const FArrayBox& /*reactfab*/, const FArrayBox& /*pressfab*/,
                       const Geometry& /*geomdata*/,
                       Real /*time*/, const Vector<BCRec>& /*bcrec*/, int /*level*/)
{
    AMREX_ASSERT(derfab.box().contains(bx));
    AMREX_ASSERT(statefab.box().contains(bx));
    AMREX_ASSERT(derfab.nComp() >= dcomp + ncomp);

    FArrayBox dummies(bx,NUM_SPECIES+1,The_Async_Arena());
    auto const& rhoY = statefab.const_array(FIRSTSPEC);
    auto const& T    = statefab.array(TEMP);
    auto       rhoD  = dummies.array(1);
    auto     lambda  = derfab.array(dcomp);
    auto         mu  = dummies.array(0);
    auto const* ltransparm = a_pelelm->trans_parms.device_trans_parm();
    amrex::ParallelFor(bx,
    [rhoY,T,rhoD,lambda,mu,ltransparm] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        getTransportCoeff(i, j, k, rhoY, T, rhoD, lambda, mu, ltransparm);
    });
}
