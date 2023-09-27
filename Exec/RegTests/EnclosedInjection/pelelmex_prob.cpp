#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", prob_parm->P_mean);
  pp.query("T_mean", prob_parm->T_mean);
  pp.query("T_jet", prob_parm->T_jet);
  pp.query("V_jet", prob_parm->V_jet);
  pp.query("jet_rad", prob_parm->jet_rad);
  pp.query("bl_thickness", prob_parm->bl_thickness);
}
