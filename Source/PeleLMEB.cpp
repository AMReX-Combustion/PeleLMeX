#include <PeleLM.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <hydro_redistribution.H>

using namespace amrex;

void PeleLM::makeEBGeometry()
{
    // TODO extend
    int max_coarsening_level = 100;
    int req_coarsening_level = 2;
    EB2::Build(geom.back(),req_coarsening_level,max_coarsening_level);
}

void PeleLM::redistributeAofS(int a_lev,
                              Real &a_dt,
                              MultiFab &a_tmpDiv, int div_comp,
                              MultiFab &a_AofS, int aofs_comp,
                              MultiFab &a_state, int state_comp,
                              int ncomp,
                              const BCRec * d_bc,
                              const Geometry &a_geom)
{
    AMREX_ASSERT(a_tmpDiv.nComp() >= div_comp+ncomp);
    AMREX_ASSERT(a_AofS.nComp() >= aofs_comp+ncomp);
    AMREX_ASSERT(a_state.nComp() >= state_comp+ncomp);

    //----------------------------------------------------------------
    const auto& ebfact = EBFactory(a_lev);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(a_tmpDiv,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

        Box const& bx = mfi.tilebox();
        auto const& flagfab = ebfact.getMultiEBCellFlagFab()[mfi];
        auto const& flag    = flagfab.const_array();
        auto const& divT_ar = a_tmpDiv.array(mfi,div_comp);
        auto const& aofs_ar = a_AofS.array(mfi,aofs_comp);

        if (flagfab.getType(bx) != FabType::covered ) {
            if (flagfab.getType(grow(bx,4)) != FabType::regular) {
                AMREX_D_TERM( auto apx = ebfact.getAreaFrac()[0]->const_array(mfi);,
                              auto apy = ebfact.getAreaFrac()[1]->const_array(mfi);,
                              auto apz = ebfact.getAreaFrac()[2]->const_array(mfi); );
                AMREX_D_TERM( Array4<Real const> fcx = ebfact.getFaceCent()[0]->const_array(mfi);,
                              Array4<Real const> fcy = ebfact.getFaceCent()[1]->const_array(mfi);,
                              Array4<Real const> fcz = ebfact.getFaceCent()[2]->const_array(mfi););
                auto const &ccc       = ebfact.getCentroid().const_array(mfi);
                auto const &vfrac_arr = ebfact.getVolFrac().const_array(mfi);

                // This is scratch space if calling StateRedistribute,
                // but is used as the weights (here set to 1) if calling
                // FluxRedistribute
                Box gbx = bx;

                if (m_adv_redist_type == "StateRedist" || 
                    m_adv_redist_type == "NewStateRedist")
                  gbx.grow(3);
                else if (m_adv_redist_type == "FluxRedist")
                  gbx.grow(2);

                FArrayBox tmpfab(gbx, ncomp);
                Elixir eli = tmpfab.elixir();
                Array4<Real> scratch = tmpfab.array(0);
                if (m_adv_redist_type == "FluxRedist")
                {    
                    amrex::ParallelFor(Box(scratch),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { scratch(i,j,k) = 1.;});   // TODO might want to test volfrac
                }
                Redistribution::Apply( bx, ncomp, aofs_ar, divT_ar,
                                       a_state.const_array(mfi, state_comp), scratch, flag,
                                       AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                       AMREX_D_DECL(fcx,fcy,fcz), ccc, d_bc,
                                       a_geom, a_dt, m_adv_redist_type );
            } else {
                // Move data to AofS for regular bx
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { aofs_ar(i,j,k,n) = divT_ar(i,j,k,n);});
            }
        }
    }
}

#endif
