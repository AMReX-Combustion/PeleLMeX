#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <pmf.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean", prob_parm->P_mean);
   pp.query("T_mean", prob_parm->T_mean);
   pp.query("T_bubble", prob_parm->T_bubble);
   pp.query("bubble_radius", prob_parm->bubble_rad);
   pp.query("bubble_y0", prob_parm->bubble_y0);
}
