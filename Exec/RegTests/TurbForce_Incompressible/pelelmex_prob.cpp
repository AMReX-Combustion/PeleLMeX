#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("turbforce");
  pp.query("urms", PeleLM::prob_parm->urms);
}

void
PeleLM::freeProbParm()
{
}
