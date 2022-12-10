#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");

   // CF params
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("T_oxidizer", PeleLM::prob_parm->T_ox);
   pp.query("T_fuel", PeleLM::prob_parm->T_fuel);
   pp.query("massflow", PeleLM::prob_parm->mdot);
   pp.query("pertmag",  PeleLM::prob_parm->pertmag);
   pp.query("jet_radius",  PeleLM::prob_parm->jetRadius);
   pp.query("inert_radius",  PeleLM::prob_parm->inertRadius);
   pp.query("inert_velocity",  PeleLM::prob_parm->inertVel);

   // ignition params
   pp.query("do_ignition",  PeleLM::prob_parm->do_ignit);
   pp.query("ignition_SphRad",  PeleLM::prob_parm->ignitSphereRad);
   pp.query("ignition_SphT",  PeleLM::prob_parm->ignitT);

   // spray params
   pp.query("mass_flow_rate",PeleLM::prob_parm->mass_flow_rate);
   pp.query("jet_vel",PeleLM::prob_parm->jet_vel);
   pp.query("jet_dx_mod",PeleLM::prob_parm->jet_dx_mod);
   pp.query("part_mean_dia",PeleLM::prob_parm->part_mean_dia);
   pp.query("part_stdev_dia",PeleLM::prob_parm->part_stdev_dia);
   pp.query("part_temp",PeleLM::prob_parm->part_temp);
}
