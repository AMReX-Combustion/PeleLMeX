#include <PeleLMeX.H>

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
    pp.get("diagnostics", diags[n], n);
    std::string diag_prefix = pele_prefix + "." + diags[n];
    ParmParse ppd(diag_prefix);
    std::string diag_type;
    ppd.get("type", diag_type);
    m_diagnostics[n] = DiagBase::create(diag_type);
    m_diagnostics[n]->init(diag_prefix, diags[n]);
    m_diagnostics[n]->addVars(m_diagVars);
  }

  // Remove duplicates from m_diagVars and check that all the variables exists
  std::sort(m_diagVars.begin(), m_diagVars.end());
  auto last = std::unique(m_diagVars.begin(), m_diagVars.end());
  m_diagVars.erase(last, m_diagVars.end());
  for (auto& v : m_diagVars) {
    bool itexists =
      derive_lst.canDerive(v) || isStateVariable(v) || isReactVariable(v);
    if (!itexists) {
      Abort("Field " + v + " is not available");
    } else {
      if (derive_lst.canDerive(v)) {
        const PeleLMDeriveRec* rec = derive_lst.get(v);
        if (rec->variableComp(v) < 0) {
          std::string errmsg = "Diagnostics can't handle derived with more "
                               "than 1 component at the moment.\n";
          errmsg += "Add the desired components individually.\n";
          Abort(errmsg);
        }
      }
    }
  }
}

void
PeleLM::updateDiagnostics()
{
  // Might need to update some internal data as the grid changes
  for (const auto& m_diagnostic : m_diagnostics) {
    if (m_diagnostic->needUpdate()) {
      m_diagnostic->prepare(
        finestLevel() + 1, Geom(0, finestLevel()), boxArray(0, finestLevel()),
        dmap, m_diagVars);
    }
  }
}

void
PeleLM::doDiagnostics()
{
  BL_PROFILE("PeleLMeX::doDiagnostics()");
  // Assemble a vector of MF containing the requested data
  Vector<std::unique_ptr<MultiFab>> diagMFVec(finestLevel() + 1);
  for (int lev{0}; lev <= finestLevel(); ++lev) {
    diagMFVec[lev] =
      std::make_unique<MultiFab>(grids[lev], dmap[lev], m_diagVars.size(), 1);
    for (int v{0}; v < m_diagVars.size(); ++v) {
      std::unique_ptr<MultiFab> mf;
      mf = derive(m_diagVars[v], m_cur_time, lev, 1);
      // If the variable is a derive component, get its index from the derive
      // multifab
      // TODO: if multiple diagVars are components of the same derive, they get
      // redundantly derived each time
      int mf_idx = 0;
      const PeleLMDeriveRec* rec = derive_lst.get(m_diagVars[v]);
      if (rec != nullptr) {
        mf_idx = rec->variableComp(m_diagVars[v]);
      }
      MultiFab::Copy(*diagMFVec[lev], *mf, mf_idx, v, 1, 1);
    }
  }

  for (const auto& m_diagnostic : m_diagnostics) {
    if (m_diagnostic->doDiag(m_cur_time, m_nstep)) {
      m_diagnostic->processDiag(
        m_nstep, m_cur_time, GetVecOfConstPtrs(diagMFVec), m_diagVars);
    }
  }
}
