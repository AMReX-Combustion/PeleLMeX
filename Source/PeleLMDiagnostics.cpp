#include <PeleLM.H>

using namespace amrex;

void
PeleLM::createDiagnostics()
{
    std::string pele_prefix = "peleLM";
    ParmParse pp(pele_prefix);
    int n_diags = 0;
    n_diags = pp.countval("diagnostics");
    Vector<std::string> diags;
    if (n_diags > 0) {
        m_diagnostics.resize(n_diags);
        diags.resize(n_diags);
    }
    for (int n = 0; n < n_diags; ++n) {
        pp.get("diagnostics",diags[n],n);
        std::string diag_prefix = pele_prefix + "." + diags[n];
        ParmParse ppd(diag_prefix);
        std::string diag_type; ppd.get("type", diag_type);
        m_diagnostics[n] = DiagBase::create(diag_type);
        m_diagnostics[n]->init(diag_prefix);
    }
}

void
PeleLM::updateDiagnostics()
{
    // Might need to update some internal data as the grid changes
}

void
PeleLM::doDiagnostics()
{
    for (int n = 0; n < m_diagnostics.size(); ++n) {
        if ( m_diagnostics[n]->doDiag(m_cur_time, m_nstep) ) {
            // Might internalize the check into processDiag()
        }
    }
}
