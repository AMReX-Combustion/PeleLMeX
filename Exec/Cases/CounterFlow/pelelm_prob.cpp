#include <PeleLM.H>
#include <AMReX_ParmParse.H>
#include <pmf.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("T_oxidizer", PeleLM::prob_parm->T_ox);
   pp.query("T_fuel", PeleLM::prob_parm->T_fuel);
   pp.query("massflow", PeleLM::prob_parm->mdot);
   pp.query("pertmag",  PeleLM::prob_parm->pertmag);
   pp.query("jet_radius",  PeleLM::prob_parm->jetRadius);
   pp.query("inert_radius",  PeleLM::prob_parm->inertRadius);
   pp.query("inert_velocity",  PeleLM::prob_parm->inertVel);

   pp.query("do_ignition",  PeleLM::prob_parm->do_ignit);
   pp.query("ignition_SphRad",  PeleLM::prob_parm->ignitSphereRad);
   pp.query("ignition_SphT",  PeleLM::prob_parm->ignitT);
   std::string pmf_datafile;
   pp.query("pmf_datafile", pmf_datafile);
   int pmf_do_average = 0;
   PMF::read_pmf(pmf_datafile, pmf_do_average);
}
