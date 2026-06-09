#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");
  pp.query("U_lid", prob_parm->U_lid);
}

void
PeleLM::freeProbParm()
{
}
