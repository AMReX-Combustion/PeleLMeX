#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <pmf.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean", prob_parm->P_mean);
   pp.query("T_mean", prob_parm->T_mean);
   pp.query("T_jet", prob_parm->T_jet);
   pp.query("V_jet", prob_parm->V_jet);
   pp.query("jet_rad", prob_parm->jet_rad);
   pp.query("bl_thickness", prob_parm->bl_thickness);

   // auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
   // amrex::ParmParse pptr("transport");
   // pp.query("const_viscosity", trans_parm.const_viscosity);
   // pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
   // pp.query("const_conductivity", trans_parm.const_conductivity);
   // pp.query("const_diffusivity", trans_parm.const_diffusivity);
   // PeleLM::trans_parms.sync_to_device();
}
