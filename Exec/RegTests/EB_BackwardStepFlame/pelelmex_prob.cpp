#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("Y_fuel", prob_parm->Y_fuel);
  pp.query("Y_oxid", prob_parm->Y_o2);
  pp.query("T_hot", prob_parm->T_hot);
  pp.query("T_wall", prob_parm->Twall);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);
  pp.query("do_ignite", prob_parm->ignite_flow);

  if (prob_parm->ignite_flow) {
    pp.query("ignition_temperature", prob_parm->ignition_temperature);
  }

  pp.query("EBinflow_T", prob_parm->EBinflow_T);
  pp.query("EBinflow_vel", prob_parm->EBinflow_vel);
  pp.query("EBinflow_Yfuel", prob_parm->EBinflow_Yfuel);

  amrex::ParmParse pplm("peleLM");
  pplm.query("EBinflow", prob_parm->do_EBinflow);
}

void
PeleLM::freeProbParm()
{
}
