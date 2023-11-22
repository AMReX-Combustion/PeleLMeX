#include <PeleLMeX.H>
#include <pelelmex_prob.H>

void
PeleLM::readProbParm()
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("standoff", PeleLM::prob_parm->standoff);
  PeleLM::pmf_data.initialize();
  amrex::Real moments[NUM_SOOT_MOMENTS + 1] = {0.0};
  if (PeleLM::do_soot_solve) {
    SootData* const sd = PeleLM::soot_model->getSootData();
    sd->initialSmallMomVals(moments);
  }
  for (int n = 0; n < NUM_SOOT_MOMENTS + 1; ++n) {
    PeleLM::prob_parm->soot_vals[n] = moments[n];
  }
}
