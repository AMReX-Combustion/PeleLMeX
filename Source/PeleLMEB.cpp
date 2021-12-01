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
    int req_coarsening_level = 0;

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

                if (m_diff_redist_type == "StateRedist" || 
                    m_diff_redist_type == "NewStateRedist")
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
    coveredState_h.resize(NVAR);
    AMREX_D_TERM(coveredState_h[0] = 0.0;,
                 coveredState_h[1] = 0.0;,
                 coveredState_h[2] = 0.0;)
    coveredState_h[DENSITY] = 1.179;
    int idO2 = O2_ID;
    int idN2 = N2_ID;
    for (int n = 0; n < NUM_SPECIES; n++ ) {
       if ( n == idO2 ) {
          coveredState_h[FIRSTSPEC+n] = 0.233*1.179;
       } else if ( n == idN2 ) {
          coveredState_h[FIRSTSPEC+n] = 0.767*1.179;
       } else {
          coveredState_h[FIRSTSPEC+n] = 0.0;
       }
    }
    coveredState_h[RHOH] = -139.7;
    coveredState_h[TEMP] = 300.0;
    coveredState_h[RHORT] = 101325.0;

    coveredState_d.resize(NVAR);
#ifdef AMREX_USE_GPU
    Gpu::htod_memcpy
#else
    std::memcpy
#endif
       (coveredState_d.data(),coveredState_h.data(), sizeof(Real)*NVAR);
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

    // TODO: do a single kernel launch there
    EB_set_covered(ldata_p->velocity,coveredState_h[0]);
    EB_set_covered(ldata_p->density,coveredState_h[DENSITY]);
    Vector<Real> covSpec(NUM_SPECIES);
    for (int n = 0; n<NUM_SPECIES; n++) {
        covSpec[n] = coveredState_h[FIRSTSPEC+n];
    }
    EB_set_covered(ldata_p->species,0,NUM_SPECIES,covSpec);
    EB_set_covered(ldata_p->rhoh,coveredState_h[RHOH]);
    EB_set_covered(ldata_p->temp,coveredState_h[TEMP]);
    EB_set_covered(ldata_p->rhoRT,coveredState_h[RHORT]);
}

void PeleLM::initialRedistribution()
{
    // Redistribute the initial solution if adv/diff scheme uses State or NewState
    if (m_adv_redist_type == "StateRedist" || 
        m_adv_redist_type == "NewStateRedist" ||
        m_diff_redist_type == "StateRedist" ||
        m_diff_redist_type == "NewStateRedist") { 

        for (int lev = 0; lev <= finest_level; ++lev) {

            // New time
            Real timeNew = getTime(lev, AmrNewTime);

            // Jugle with Old/New states: fillPatch old and redistribute
            // from Old->New to end up with proper New state
            auto ldataOld_p = getLevelDataPtr(lev,AmrOldTime);
            auto ldataNew_p = getLevelDataPtr(lev,AmrNewTime);

            auto const& fact = EBFactory(lev);

            // Velocities
            EB_set_covered(ldataNew_p->velocity,0.0);
            ldataNew_p->velocity.FillBoundary(geom[lev].periodicity());
            MultiFab::Copy(ldataOld_p->velocity, ldataNew_p->velocity, 0, 0, AMREX_SPACEDIM, m_nGrowState);
            fillpatch_velocity(lev, timeNew, ldataOld_p->velocity, m_nGrowState);

            if (!m_incompressible) {

                // Mass
                EB_set_covered(ldataNew_p->density,0.0);
                ldataNew_p->density.FillBoundary(geom[lev].periodicity());
                MultiFab::Copy(ldataOld_p->density, ldataNew_p->density, 0, 0, 1, m_nGrowState);
                fillpatch_density(lev, timeNew, ldataOld_p->density, m_nGrowState);

                EB_set_covered(ldataNew_p->species,0.0);
                ldataNew_p->species.FillBoundary(geom[lev].periodicity());
                MultiFab::Copy(ldataOld_p->species, ldataNew_p->species, 0, 0, NUM_SPECIES, m_nGrowState);
                fillpatch_species(lev, timeNew, ldataOld_p->species, m_nGrowState);

                // Energy
                EB_set_covered(ldataNew_p->rhoh,0.0);
                ldataNew_p->rhoh.FillBoundary(geom[lev].periodicity());
                MultiFab::Copy(ldataOld_p->rhoh, ldataNew_p->rhoh, 0, 0, 1, m_nGrowState);
                EB_set_covered(ldataNew_p->temp,0.0);
                ldataNew_p->temp.FillBoundary(geom[lev].periodicity());
                MultiFab::Copy(ldataOld_p->temp, ldataNew_p->temp, 0, 0, 1, m_nGrowState);
                fillpatch_energy(lev, timeNew, ldataOld_p->rhoh, ldataOld_p->temp, m_nGrowState);
            }

            for (MFIter mfi(ldataNew_p->velocity,TilingIfNotGPU()); mfi.isValid(); ++mfi) {   

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

                    auto bcRecVel = fetchBCRecArray(VELX,AMREX_SPACEDIM); 
                    auto bcRecVel_d = convertToDeviceVector(bcRecVel);
                    Redistribution::ApplyToInitialData( bx, AMREX_SPACEDIM,
                                                        ldataNew_p->velocity.array(mfi), ldataOld_p->velocity.array(mfi),
                                                        flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                        AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                        bcRecVel_d.dataPtr(), geom[lev], m_adv_redist_type);
                    if (!m_incompressible) {
                        auto bcRec = fetchBCRecArray(FIRSTSPEC,NUM_SPECIES); 
                        auto bcRec_d = convertToDeviceVector(bcRec);
                        Redistribution::ApplyToInitialData( bx, 1,
                                                            ldataNew_p->density.array(mfi), ldataOld_p->density.array(mfi),
                                                            flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                            AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                            bcRec_d.dataPtr(), geom[lev], m_adv_redist_type);
                        Redistribution::ApplyToInitialData( bx, NUM_SPECIES,
                                                            ldataNew_p->species.array(mfi), ldataOld_p->species.array(mfi),
                                                            flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                            AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                            bcRec_d.dataPtr(), geom[lev], m_adv_redist_type);
                        auto bcRecRhoH = fetchBCRecArray(RHOH,1); 
                        auto bcRecRhoH_d = convertToDeviceVector(bcRecRhoH);
                        Redistribution::ApplyToInitialData( bx, 1,
                                                            ldataNew_p->rhoh.array(mfi), ldataOld_p->rhoh.array(mfi),
                                                            flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                            AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                            bcRecRhoH_d.dataPtr(), geom[lev], m_adv_redist_type);
                        auto bcRecTemp = fetchBCRecArray(TEMP,1); 
                        auto bcRecTemp_d = convertToDeviceVector(bcRecTemp);
                        Redistribution::ApplyToInitialData( bx, 1,
                                                            ldataNew_p->temp.array(mfi), ldataOld_p->temp.array(mfi),
                                                            flag, AMREX_D_DECL(apx, apy, apz), vfrac,
                                                            AMREX_D_DECL(fcx, fcy, fcz), ccc,
                                                            bcRecTemp_d.dataPtr(), geom[lev], m_adv_redist_type);
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
#endif
