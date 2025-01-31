#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  // CF params
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("T_oxidizer", PeleLM::prob_parm->T_ox);
  pp.query("T_fuel", PeleLM::prob_parm->T_fuel);
  pp.query("T_inert", PeleLM::prob_parm->T_inert);
  pp.query("Y_O2_ox", PeleLM::prob_parm->Y_O2_ox);
  pp.query("Y_N2_ox", PeleLM::prob_parm->Y_N2_ox);
  pp.query("Y_N2_fuel", PeleLM::prob_parm->Y_N2_fuel);
  pp.query("Y_fuel", PeleLM::prob_parm->Y_fuel);
  pp.query("Y_H2O_ign", PeleLM::prob_parm->Y_H2O_ign);
  pp.query("Y_CO2_ign", PeleLM::prob_parm->Y_CO2_ign);
  pp.query("massflow_ox", PeleLM::prob_parm->mdot_ox);
  pp.query("massflow_fuel", PeleLM::prob_parm->mdot_fuel);
  pp.query("pertmag", PeleLM::prob_parm->pertmag);
  pp.query("jet_radius", PeleLM::prob_parm->jetRadius);
  pp.query("inert_radius", PeleLM::prob_parm->inertRadius);
  pp.query("inert_velocity_ox", PeleLM::prob_parm->inertVel_ox);
  pp.query("inert_velocity_fuel", PeleLM::prob_parm->inertVel_fuel);

  // ignition params
  pp.query("do_ignition", PeleLM::prob_parm->do_ignit);
  pp.query("ignition_SphRad", PeleLM::prob_parm->ignitSphereRad);
  pp.query("ignition_SphT", PeleLM::prob_parm->ignitT);
}
