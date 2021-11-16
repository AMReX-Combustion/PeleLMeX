#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("T_mean", PeleLM::prob_parm->T_mean);
   pp.query("rvort", PeleLM::prob_parm->rvort);
   pp.query("xvort", PeleLM::prob_parm->xvort);
   pp.query("yvort", PeleLM::prob_parm->yvort);
   pp.query("forcevort", PeleLM::prob_parm->forcevort);
   pp.query("centx", PeleLM::prob_parm->centx);
   pp.query("centy", PeleLM::prob_parm->centy);
   pp.query("r_circ", PeleLM::prob_parm->r_circ);
   pp.query("r_hole", PeleLM::prob_parm->r_hole);
   pp.query("nholes", PeleLM::prob_parm->nholes);
   pp.query("cone_angle", PeleLM::prob_parm->cone_angle);
   pp.query("T_jet", PeleLM::prob_parm->T_jet);
   pp.query("vel_jet", PeleLM::prob_parm->vel_jet);

   /*
   if (!m_incompressible) {
      auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
      amrex::ParmParse pptr("transport");
      pp.query("const_viscosity", trans_parm.const_viscosity);
      pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
      pp.query("const_conductivity", trans_parm.const_conductivity);
      pp.query("const_diffusivity", trans_parm.const_diffusivity);
      PeleLM::trans_parms.sync_to_device();
   }
   */
}
