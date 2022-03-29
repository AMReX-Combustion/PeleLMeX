#include "DiagFramePlane.H"

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
                        const amrex::Vector<amrex::BoxArray> &a_grids)
{
}

void
DiagFramePlane::processDiag(const amrex::Real &a_time,
                            const amrex::Vector<const amrex::MultiFab*> &a_state)
{
}
