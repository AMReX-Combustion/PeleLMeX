#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  std::string type;
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  pp.query("pertmag", PeleLM::prob_parm->pertmag);
  pp.query("Vin", PeleLM::prob_parm->Vin);
  pp.query("Vcoflow", PeleLM::prob_parm->Vcoflow);
  pp.query("T_lips", PeleLM::prob_parm->T_lips);
  pp.query("r_int", PeleLM::prob_parm->r_int);
  pp.query("r_ext", PeleLM::prob_parm->r_ext);
  PeleLM::prob_parm->phi = -1;
  pp.query("equivalence_ratio", PeleLM::prob_parm->phi);
  if (PeleLM::prob_parm->phi == 1) {
    amrex::Print() << " eq_ratio is stoichiometric \n";
  } else if (PeleLM::prob_parm->phi == 2) {
    amrex::Print() << " eq_ratio is lean \n";
  } else if (PeleLM::prob_parm->phi == 3) {
    amrex::Print() << " eq_ratio is rich \n";
  } else {
    amrex::Abort("Not a valid equivalence_ratio");
  }

#ifdef PELE_USE_PLASMA
  pp.query("electrode_radius", PeleLM::prob_parm->electrode_radius);
  pp.query("electrode_width", PeleLM::prob_parm->electrode_width);
  pp.query("electrode_phiV", PeleLM::prob_parm->electrode_phiV);
  pp.query("burner_phiV", PeleLM::prob_parm->burner_phiV);
  pp.query("burner_Rext", PeleLM::prob_parm->burner_Rext);
#endif

  PeleLM::pmf_data.initialize();
}
