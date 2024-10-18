#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("meanFlowDir", prob_parm->meanFlowDir);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);
  pp.query("ode_IC", prob_parm->ode_IC);
  pp.query("ode_xy_lo", prob_parm->ode_xy_lo);
  pp.query("ode_length", prob_parm->ode_length);
  pp.query("ode_height", prob_parm->ode_height);
  pp.query("ode_srcstrength", prob_parm->ode_srcstrength);
}
