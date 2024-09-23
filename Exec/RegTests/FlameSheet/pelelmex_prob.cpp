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
  pp.query("pertlength", PeleLM::prob_parm->pertlength);

  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();
  PeleLM::pmf_data.initialize();
}
