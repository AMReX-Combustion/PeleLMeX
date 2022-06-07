#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("T_mean", prob_parm->T_mean);
   pp.query("P_mean", prob_parm->P_mean);
   pp.query("meanFlowDir", prob_parm->meanFlowDir);
   pp.query("meanFlowMag", prob_parm->meanFlowMag);
}
