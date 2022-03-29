#include "DiagBase.H"

bool
DiagBase::doDiag(const amrex::Real &a_time,
                 int a_nstep)
{
    bool willDo = false;
    if ( m_freq > 0 && (a_nstep % m_freq == 0) ) {
        willDo = true;
    }

    // TODO: using a_time

    return willDo;
}
