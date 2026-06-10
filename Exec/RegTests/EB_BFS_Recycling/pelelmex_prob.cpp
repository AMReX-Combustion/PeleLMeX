#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");
  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);
  pp.query("turb_seed_amplitude", prob_parm->turb_seed_amplitude);
}

void
PeleLM::freeProbParm()
{
}
