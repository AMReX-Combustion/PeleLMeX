#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T_mean", PeleLM::prob_parm->T_mean);
  pp.query("rvort", PeleLM::prob_parm->rvort);
  pp.query("xvort", PeleLM::prob_parm->xvort);
  pp.query("yvort", PeleLM::prob_parm->yvort);
  pp.query("forcevort", PeleLM::prob_parm->forcevort);

  if (!m_incompressible) {
    auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
    amrex::ParmParse pptr("transport");
    pp.query("const_viscosity", trans_parm.const_viscosity);
    pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
    pp.query("const_conductivity", trans_parm.const_conductivity);
    pp.query("const_diffusivity", trans_parm.const_diffusivity);
    PeleLM::trans_parms.sync_to_device();
  }
}
