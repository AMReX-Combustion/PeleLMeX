#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean", prob_parm->P_mean);
   pp.query("T_mean", prob_parm->T_mean);
   pp.query("T_bubble", prob_parm->T_bubble);
   pp.query("bubble_radius", prob_parm->bubble_rad);
   pp.query("bubble_y0", prob_parm->bubble_y0);
   pp.query("use_symmetry", prob_parm->is_sym);
   pp.query("use_mix_bubble", prob_parm->bubble_is_mix);

   auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
   amrex::ParmParse pptr("transport");
   pp.query("const_viscosity", trans_parm.const_viscosity);
   pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
   pp.query("const_conductivity", trans_parm.const_conductivity);
   pp.query("const_diffusivity", trans_parm.const_diffusivity);
   PeleLM::trans_parms.sync_to_device();
}
