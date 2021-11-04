#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   std::string type;
   pp.query("P_mean",   PeleLM::prob_parm->P_mean);
   pp.query("standoff", PeleLM::prob_parm->standoff);
   pp.query("pertmag",  PeleLM::prob_parm->pertmag);

   PeleLM::pmf_data.initialize();

   auto& trans_parm = PeleLM::trans_parms.host_trans_parm();
   amrex::ParmParse pptr("transport");
   pp.query("const_viscosity", trans_parm.const_viscosity);
   pp.query("const_bulk_viscosity", trans_parm.const_bulk_viscosity);
   pp.query("const_conductivity", trans_parm.const_conductivity);
   pp.query("const_diffusivity", trans_parm.const_diffusivity);
   PeleLM::trans_parms.sync_to_device();
}
