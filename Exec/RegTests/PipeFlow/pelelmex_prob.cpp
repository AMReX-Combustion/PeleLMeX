#include <PeleLMeX.H>
#include <AMReX_ParmParse.H>

void
PeleLM::readProbParm() // NOLINT(readability-make-member-function-const)
{
  amrex::ParmParse pp("prob");

  pp.query("T_mean", prob_parm->T_mean);
  pp.query("P_mean", prob_parm->P_mean);
  pp.query("meanFlowDir", prob_parm->meanFlowDir);
  pp.query("meanFlowMag", prob_parm->meanFlowMag);
  pp.query("perturbMag", prob_parm->perturbMag);
  pp.query("problem_type", prob_parm->flowType);
  AMREX_ALWAYS_ASSERT(prob_parm->flowType == 1 || prob_parm->flowType == 2);

  // If mesh mapping is active, pick up the ConstantMap scaling factors
  // so the IC can be evaluated in physical coordinates.  Falls back to
  // fac = 1 when mapping is off.
  {
    amrex::ParmParse ppcm("ConstantMap");
    amrex::Vector<amrex::Real> fac(AMREX_SPACEDIM, 1.0);
    ppcm.queryarr("scaling_factor", fac, 0, AMREX_SPACEDIM);
    prob_parm->fac_x = fac[0];
#if AMREX_SPACEDIM >= 2
    prob_parm->fac_y = fac[1];
#endif
#if AMREX_SPACEDIM >= 3
    prob_parm->fac_z = fac[2];
#endif
  }
}

void
PeleLM::freeProbParm()
{
}
