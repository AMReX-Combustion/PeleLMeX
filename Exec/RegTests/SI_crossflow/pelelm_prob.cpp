#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
    amrex::ParmParse pp("prob");
   
    pp.query("T_mean", prob_parm->T_mean);
    pp.query("P_mean", prob_parm->P_mean);
    pp.query("flowDir", prob_parm->meanFlowDir);
    pp.query("flowMag", prob_parm->meanFlowMag);

    pp.query("V_in"  , prob_parm->V_in);

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
