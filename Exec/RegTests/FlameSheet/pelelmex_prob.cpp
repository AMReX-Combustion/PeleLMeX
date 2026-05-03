#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  pp.query("pertmag", PeleLM::prob_parm->pertmag);
  pp.query("pertlength", PeleLM::prob_parm->pertlength);
  pp.query("ref_frame", PeleLM::prob_parm->ref_frame);
  if (PeleLM::prob_parm->ref_frame == PmfReferenceFrame::fixed_velocity) {
    pp.get("moving_frame_velocity", PeleLM::prob_parm->moving_frame_velocity);
  }

  PeleLM::prob_parm->eosparm = PeleLM::eos_parms.device_parm();
  PeleLM::pmf_data.initialize();
}

void
PeleLM::freeProbParm()
{
  PeleLM::pmf_data.deallocate();
}
