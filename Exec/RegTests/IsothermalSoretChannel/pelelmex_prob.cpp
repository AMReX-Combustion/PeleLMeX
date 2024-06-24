#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);

  // TODO: somewhat hard coded bath, fuel and oxid IDs
  // should exist somewhere in PeleLM.
  PeleLM::prob_parm->bathID = N2_ID;
  PeleLM::prob_parm->fuelID = H2_ID;
  PeleLM::prob_parm->oxidID = O2_ID;
}
