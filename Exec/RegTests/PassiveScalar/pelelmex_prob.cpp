#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");
  pp.query("P_mean", prob_parm->P_mean);
  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();
}
