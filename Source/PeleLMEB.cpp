#include <PeleLM.H>
#include <PeleLMUtils.H>
#include <EBUserDefined.H>

#ifdef AMREX_USE_EB
#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <hydro_redistribution.H>
#include <AMReX_EBMFInterpolater.H>

using namespace amrex;

void PeleLM::makeEBGeometry()
{
    // TODO extend
    int max_coarsening_level = 100;
    int req_coarsening_level = geom.size()-1;

    // Read the geometry type and act accordingly
    ParmParse ppeb2("eb2");
    std::string geom_type;
    ppeb2.get("geom_type", geom_type);

    if (geom_type == "UserDefined") {
        EBUserDefined(geom.back(), req_coarsening_level, max_coarsening_level);
    } else {
        // If geom_type is not an AMReX recognized type, it'll crash.
        EB2::Build(geom.back(), req_coarsening_level, max_coarsening_level);
    }
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

                if (m_adv_redist_type == "StateRedist")
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

void PeleLM::redistributeDiff(int a_lev,
                              const Real &a_dt,
                              MultiFab &a_tmpDiv, int div_comp,
                              MultiFab &a_diff, int diff_comp,
                              const MultiFab &a_state, int state_comp,
                              int ncomp,
                              const BCRec * d_bc,
                              const Geometry &a_geom)
{
    AMREX_ASSERT(a_tmpDiv.nComp() >= div_comp+ncomp);
    AMREX_ASSERT(a_diff.nComp() >= diff_comp+ncomp);
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
        auto const& diff_ar = a_diff.array(mfi,diff_comp);

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

                if (m_diff_redist_type == "StateRedist")
                  gbx.grow(3);
                else if (m_diff_redist_type == "FluxRedist")
                  gbx.grow(2);

                FArrayBox tmpfab(gbx, ncomp);
                Elixir eli = tmpfab.elixir();
                Array4<Real> scratch = tmpfab.array(0);
                if (m_diff_redist_type == "FluxRedist")
                {
                    amrex::ParallelFor(Box(scratch),
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    { scratch(i,j,k) = 1.;});   // TODO might want to test volfrac
                }
                Redistribution::Apply( bx, ncomp, diff_ar, divT_ar,
                                       a_state.const_array(mfi, state_comp), scratch, flag,
                                       AMREX_D_DECL(apx,apy,apz), vfrac_arr,
                                       AMREX_D_DECL(fcx,fcy,fcz), ccc, d_bc,
                                       a_geom, a_dt, m_diff_redist_type );
            } else {
                // Move data to AofS for regular bx
                amrex::ParallelFor(bx, ncomp, [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                { diff_ar(i,j,k,n) = divT_ar(i,j,k,n);});
            }
        }
    }
}

void PeleLM::initCoveredState()
{
    // TODO use typical values
    if ( m_incompressible ) {
        coveredState_h.resize(AMREX_SPACEDIM);
        AMREX_D_TERM(coveredState_h[0] = typical_values[VELX+0];,
                     coveredState_h[1] = typical_values[VELX+1];,
                     coveredState_h[2] = typical_values[VELX+2];)
        coveredState_d.resize(AMREX_SPACEDIM);
#ifdef AMREX_USE_GPU
        Gpu::htod_memcpy
#else
        std::memcpy
#endif
          (coveredState_d.data(),coveredState_h.data(), sizeof(Real)*AMREX_SPACEDIM);
    } else {
        coveredState_h.resize(NVAR);
        AMREX_D_TERM(coveredState_h[0] = typical_values[VELX+0];,
                     coveredState_h[1] = typical_values[VELX+1];,
                     coveredState_h[2] = typical_values[VELX+2];)
        coveredState_h[DENSITY] = typical_values[DENSITY];
        for (int n = 0; n < NUM_SPECIES; n++ ) {
           coveredState_h[FIRSTSPEC+n] = typical_values[FIRSTSPEC+n];
        }
        coveredState_h[RHOH] = typical_values[RHOH];
        coveredState_h[TEMP] = typical_values[TEMP]-10.0;
        coveredState_h[RHORT] = typical_values[RHORT];

        coveredState_d.resize(NVAR);
#ifdef AMREX_USE_GPU
        Gpu::htod_memcpy
#else
        std::memcpy
#endif
          (coveredState_d.data(),coveredState_h.data(), sizeof(Real)*NVAR);
    }
}

void PeleLM::setCoveredState(const TimeStamp &a_time)
{
    for (int lev = 0; lev <= finest_level; lev++) {
        setCoveredState(lev,a_time);
    }
}

void PeleLM::setCoveredState(int lev, const TimeStamp &a_time)
{
    AMREX_ASSERT(a_time == AmrOldTime || a_time == AmrNewTime);

    auto ldata_p = getLevelDataPtr(lev,a_time);

    if ( m_incompressible ) {
        EB_set_covered(ldata_p->state,0,AMREX_SPACEDIM,coveredState_h);
    } else {
        EB_set_covered(ldata_p->state,0,NVAR,coveredState_h);
    }
}

void PeleLM::initialRedistribution()
{
    // Redistribute the initial solution if adv/diff scheme uses State or NewState
    if (m_adv_redist_type == "StateRedist" ||
        m_diff_redist_type == "StateRedist") {

        for (int lev = 0; lev <= finest_level; ++lev) {

            // New time
            Real timeNew = getTime(lev, AmrNewTime);

            // Jungle with Old/New states: fillPatch old and redistribute
            // from Old->New to end up with proper New state
            auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
            auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

            auto const& fact = EBFactory(lev);

            // State
            if ( m_incompressible ) {
               Vector<Real> stateCovered(AMREX_SPACEDIM,0.0);
               EB_set_covered(ldataNew_p->state,0,AMREX_SPACEDIM,stateCovered);
               ldataNew_p->state.FillBoundary(geom[lev].periodicity());
               MultiFab::Copy(ldataOld_p->state, ldataNew_p->state, 0, 0, AMREX_SPACEDIM, m_nGrowState);
            } else {
               Vector<Real> stateCovered(NVAR,0.0);
               EB_set_covered(ldataNew_p->state,0,NVAR,stateCovered);
               ldataNew_p->state.FillBoundary(geom[lev].periodicity());
               MultiFab::Copy(ldataOld_p->state, ldataNew_p->state, 0, 0, NVAR, m_nGrowState);
            }
            fillpatch_state(lev, timeNew, ldataOld_p->state, m_nGrowState);

            for (MFIter mfi(ldataNew_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi) {

                const Box& bx = mfi.validbox();

                EBCellFlagFab const& flagfab = fact.getMultiEBCellFlagFab()[mfi];
                Array4<EBCellFlag const> const& flag = flagfab.const_array();

                if ( (flagfab.getType(bx)                != FabType::covered) &&
                     (flagfab.getType(amrex::grow(bx,4)) != FabType::regular) )
                {
                    Array4<Real const> AMREX_D_DECL(fcx, fcy, fcz), ccc, vfrac, AMREX_D_DECL(apx, apy, apz);
                    AMREX_D_TERM(fcx = fact.getFaceCent()[0]->const_array(mfi);,
                                 fcy = fact.getFaceCent()[1]->const_array(mfi);,
                                 fcz = fact.getFaceCent()[2]->const_array(mfi););
                    ccc   = fact.getCentroid().const_array(mfi);
                    AMREX_D_TERM(apx = fact.getAreaFrac()[0]->const_array(mfi);,
                                 apy = fact.getAreaFrac()[1]->const_array(mfi);,
                                 apz = fact.getAreaFrac()[2]->const_array(mfi););
                    vfrac = fact.getVolFrac().const_array(mfi);

                    if ( m_incompressible ) {
                        auto bcRec = fetchBCRecArray(0,AMREX_SPACEDIM);
                        auto bcRec_d = convertToDeviceVector(bcRec);
                        Redistribution::ApplyToInitialData( bx, AMREX_SPACEDIM,
                                                            ldataNew_p->state.array(mfi,0), ldataOld_p->state.array(mfi,0),
                                                            flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                            AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                            bcRec_d.dataPtr(), geom[lev], m_adv_redist_type);
                    } else {
                        auto bcRec = fetchBCRecArray(0,NVAR);
                        auto bcRec_d = convertToDeviceVector(bcRec);
                        Redistribution::ApplyToInitialData( bx, NVAR,
                                                            ldataNew_p->state.array(mfi,0), ldataOld_p->state.array(mfi,0),
                                                            flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                            AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                            bcRec_d.dataPtr(), geom[lev], m_adv_redist_type);
                    }
                }
            }
        }
    }

    // Initialize covered state
    setCoveredState(AmrNewTime);
}

void PeleLM::getEBDistance(int a_lev,
                           MultiFab &a_signDistLev) {


    if (a_lev == 0) {
        MultiFab::Copy(a_signDistLev,*m_signedDist0,0,0,1,0);
        return;
    }

    // A pair of MF to hold crse & fine dist data
    Array<std::unique_ptr<MultiFab>,2> MFpair;

    // Interpolate on successive levels up to a_lev
    for (int lev = 1; lev <= a_lev; ++lev) {

        // Use MF EB interp
        auto& interpolater  = eb_mf_lincc_interp;

        // Get signDist on coarsen fineBA
        BoxArray coarsenBA(grids[lev].size());
        for (int j = 0, N = coarsenBA.size(); j < N; ++j)
        {
            coarsenBA.set(j,interpolater.CoarseBox(grids[lev][j], refRatio(lev-1)));
        }
        MultiFab coarsenSignDist(coarsenBA,dmap[lev],1,0);
        coarsenSignDist.setVal(0.0);
        MultiFab *crseSignDist = (lev == 1) ? m_signedDist0.get()
                                            : MFpair[0].get();
        coarsenSignDist.ParallelCopy(*crseSignDist,0,0,1);

        // Interpolate on current lev
        MultiFab *currentSignDist;
        if ( lev < a_lev ) {
            MFpair[1].reset(new MultiFab(grids[lev],dmap[lev],1,0,MFInfo(),EBFactory(lev)));
        }
        currentSignDist = (lev == a_lev) ? &a_signDistLev
                                         : MFpair[1].get();

        interpolater.interp(coarsenSignDist, 0,
                            *currentSignDist, 0,
                            1, IntVect(0),
                            Geom(lev-1), Geom(lev),
                            Geom(lev).Domain(), refRatio(lev-1),
                            {m_bcrec_force},0);

        // Swap MFpair
        if (lev < a_lev ) {
            swap(MFpair[0],MFpair[1]);
        }
    }
}

void
PeleLM::correct_vel_small_cells (Vector<MultiFab*> const& a_vel,
                                 Vector<Array<MultiFab const*,AMREX_SPACEDIM> > const& a_umac)
{
    BL_PROFILE("PeleLM::correct_vel_small_cells");

    for (int lev = 0; lev <= finest_level; lev++)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*a_vel[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
       {
          // Tilebox
          const Box bx = mfi.tilebox();

          EBCellFlagFab const& flags = EBFactory(lev).getMultiEBCellFlagFab()[mfi];

          // Face-centered velocity components
          AMREX_D_TERM(const auto& umac_fab = (a_umac[lev][0])->const_array(mfi);,
                       const auto& vmac_fab = (a_umac[lev][1])->const_array(mfi);,
                       const auto& wmac_fab = (a_umac[lev][2])->const_array(mfi););

          if (flags.getType(amrex::grow(bx,0)) == FabType::covered )
          {
            // do nothing
          }

          // No cut cells in this FAB
          else if (flags.getType(amrex::grow(bx,1)) == FabType::regular )
          {
            // do nothing
          }

          // Cut cells in this FAB
          else
          {
             // Face-centered areas
             AMREX_D_TERM(const auto& apx_fab   = EBFactory(lev).getAreaFrac()[0]->const_array(mfi);,
                          const auto& apy_fab   = EBFactory(lev).getAreaFrac()[1]->const_array(mfi);,
                          const auto& apz_fab   = EBFactory(lev).getAreaFrac()[2]->const_array(mfi););

             const auto& vfrac_fab = EBFactory(lev).getVolFrac().const_array(mfi);

             const auto& ccvel_fab = a_vel[lev]->array(mfi);

             // This FAB has cut cells -- we define the centroid value in terms of the MAC velocities onfaces
             amrex::ParallelFor(bx,
               [vfrac_fab,AMREX_D_DECL(apx_fab,apy_fab,apz_fab),ccvel_fab,AMREX_D_DECL(umac_fab,vmac_fab,wmac_fab)]
               AMREX_GPU_DEVICE (int i, int j, int k) noexcept
             {
                 if (vfrac_fab(i,j,k) > 0.0 && vfrac_fab(i,j,k) < 5.e-3)
                 {
                    AMREX_D_TERM(Real u_avg = (apx_fab(i,j,k) * umac_fab(i,j,k) + apx_fab(i+1,j,k) * umac_fab(i+1,j,k))
                                            / (apx_fab(i,j,k) + apx_fab(i+1,j,k));,
                                 Real v_avg = (apy_fab(i,j,k) * vmac_fab(i,j,k) + apy_fab(i,j+1,k) * vmac_fab(i,j+1,k))
                                            / (apy_fab(i,j,k) + apy_fab(i,j+1,k));,
                                 Real w_avg = (apz_fab(i,j,k) * wmac_fab(i,j,k) + apz_fab(i,j,k+1) * wmac_fab(i,j,k+1))
                                            / (apz_fab(i,j,k) + apz_fab(i,j,k+1)););

                    AMREX_D_TERM(ccvel_fab(i,j,k,0) = u_avg;,
                                 ccvel_fab(i,j,k,1) = v_avg;,
                                 ccvel_fab(i,j,k,2) = w_avg;);

                 }
             });
          }
       }
    }
}

#endif
