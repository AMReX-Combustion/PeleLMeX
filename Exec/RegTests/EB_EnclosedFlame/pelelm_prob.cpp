#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   std::string type;
   pp.query("P_mean",   PeleLM::prob_parm->P_mean);
   pp.query("standoff", PeleLM::prob_parm->standoff);
   pp.query("pertmag",  PeleLM::prob_parm->pertmag);
#if (AMREX_SPACEDIM == 3)
   pp.query("pertLz",  PeleLM::prob_parm->pertLz);
#endif

   PeleLM::pmf_data.initialize();
}
