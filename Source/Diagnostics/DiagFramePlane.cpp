#include "DiagFramePlane.H"
#include "AMReX_VisMF.H"
#include "AMReX_PlotFileUtil.H"

void
DiagFramePlane::init(const std::string &a_prefix)
{
    amrex::ParmParse pp(a_prefix);
    // Plane normal
    pp.get("normal", m_normal);
    AMREX_ASSERT(m_normal>=0 && m_normal<AMREX_SPACEDIM);

    // Plane center
    amrex::Vector<amrex::Real> center;
    pp.getarr("center",center,0,pp.countval("center"));
    if (center.size() == AMREX_SPACEDIM) {
        for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
            m_center[idim] = center[idim];
        }
    } else if (center.size() == 1) {
        m_center[m_normal] = center[0];
    }
    
    pp.query("int", m_freq);
    pp.query("per", m_per);
    AMREX_ASSERT(m_freq>0 || m_per>0.0);
}

void 
DiagFramePlane::prepare(int a_nlevels,
                        const amrex::Vector<amrex::Geometry> &a_geoms,
                        const amrex::Vector<amrex::BoxArray> &a_grids,
                        const amrex::Vector<amrex::DistributionMapping> &a_dmap)
{
    // Store the level0 geometry
    m_geomLev0 = a_geoms[0];

    // Resize internal vectors
    m_intwgt.resize(a_nlevels);
    m_k0.resize(a_nlevels);

    // On each level, find the k0 where the plane lays
    // and the weight of the directionnal interpolation
    for (int lev = 0; lev < a_nlevels; lev++) {
        auto const& domain = a_geoms[lev].Domain();  
        const amrex::Real* dx     = a_geoms[lev].CellSize();
        const amrex::Real* problo = a_geoms[lev].ProbLo();
        // How many dx away from the lowest cell-center ?
        amrex::Real dist = (m_center[m_normal] - (problo[m_normal]+0.5*dx[m_normal]))/dx[m_normal];
        int k0 = static_cast<int>(std::round(dist));
        dist -= static_cast<amrex::Real>(k0);
        m_k0[lev] = k0;
        // Quadratic interp. weights on k0-1, k0, k0+1
        m_intwgt[lev][0] = 0.5 * (dist - 1.0) * (dist - 2.0);
        m_intwgt[lev][1] = dist * (2.0 - dist);
        m_intwgt[lev][2] = 0.5 * dist * (dist - 1.0);;
    }

    // Assemble the 2D slice boxArray
    m_sliceBA.resize(a_nlevels);
    m_sliceDM.resize(a_nlevels);
    m_dmConvert.resize(a_nlevels);
    for (int lev = 0; lev < a_nlevels; lev++) {
        amrex::Vector<int> pmap;
        amrex::BoxList bl(a_grids[lev].ixType());
        bl.reserve(a_grids[lev].size());
        for (int i = 0; i < a_grids[lev].size(); ++i) {
            auto cBox = a_grids[lev][i];
            amrex::IntVect ploc(AMREX_D_DECL(cBox.smallEnd(0),cBox.smallEnd(1),cBox.smallEnd(2)));
            ploc[m_normal] = m_k0[lev];
            if (cBox.contains(ploc)) {
                cBox.setRange(m_normal,0,1);
                bl.push_back(cBox);
                pmap.push_back(a_dmap[lev][i]);
                m_dmConvert[lev].push_back(i);
            }
        }
        m_sliceBA[lev].define(bl);
        m_sliceDM[lev].define(pmap);
    }
}

void
DiagFramePlane::processDiag(const amrex::Real &a_time,
                            const amrex::Vector<const amrex::MultiFab*> &a_state)
{
    amrex::Vector<amrex::MultiFab> planeData(a_state.size());
    for (int lev = 0; lev < a_state.size(); ++lev) {
        planeData[lev].define(m_sliceBA[lev], m_sliceDM[lev], a_state[0]->nComp(), 0);
        int p0 = m_k0[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(planeData[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto &bx = mfi.tilebox();
            const int state_idx = m_dmConvert[lev][mfi.index()];
            auto const& state = a_state[lev]->const_array(state_idx,0);
            auto const& plane = planeData[lev].array(mfi);
            if (m_normal == 0) {
                amrex::ParallelFor(bx, a_state[0]->nComp(), [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    plane(i,j,k,n) = m_intwgt[lev][0] * state(p0-1,j,k,n) +
                                     m_intwgt[lev][1] * state(p0  ,j,k,n) +
                                     m_intwgt[lev][2] * state(p0+1,j,k,n);
                });
            } else if (m_normal == 1) {
                amrex::ParallelFor(bx, a_state[0]->nComp(), [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    plane(i,j,k,n) = m_intwgt[lev][0] * state(i,p0-1,k,n) +
                                     m_intwgt[lev][1] * state(i,p0  ,k,n) +
                                     m_intwgt[lev][2] * state(i,p0+1,k,n);
                });
            } else if (m_normal == 2) {
                amrex::ParallelFor(bx, a_state[0]->nComp(), [=]
                AMREX_GPU_DEVICE (int i, int j, int k, int n) noexcept
                {
                    plane(i,j,k,n) = m_intwgt[lev][0] * state(i,j,p0-1,n) +
                                     m_intwgt[lev][1] * state(i,j,p0  ,n) +
                                     m_intwgt[lev][2] * state(i,j,p0+1,n);
                });
            }
        }
    }

    amrex::Vector<amrex::Geometry> pltGeoms(a_state.size());
    pltGeoms[0] = m_geomLev0;
    amrex::IntVect rref(AMREX_D_DECL(2,2,2));
    rref[m_normal] = 1;
    for (int lev = 1; lev < a_state.size(); ++lev) {
        pltGeoms[lev] = amrex::refine(pltGeoms[lev-1], rref);
    }
    amrex::WriteMLMF("test_"+std::to_string(a_time),GetVecOfConstPtrs(planeData),pltGeoms);
    
    //amrex::VisMF::Write(planeData[lev],"planeDataLev"+std::to_string(lev)+std::to_string(a_time));
    //amrex::Abort();
}
