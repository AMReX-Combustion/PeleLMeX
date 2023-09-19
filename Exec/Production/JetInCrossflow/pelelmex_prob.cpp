#include <PeleLMeX.H>

void
PeleLM::readProbParm()
{
  amrex::ParmParse pp("prob");

  std::string type;
  pp.query("P_mean", PeleLM::prob_parm->P_mean);
  pp.query("V_in", PeleLM::prob_parm->V_in);
  pp.query("jet_rad", PeleLM::prob_parm->jet_rad);
  pp.query("jet_temp", PeleLM::prob_parm->jet_temp);
  pp.query("global_eq_ratio", PeleLM::prob_parm->global_eq_ratio);
  pp.query("ox_temp", PeleLM::prob_parm->ox_temp);
  pp.query("X_O2", PeleLM::prob_parm->X_O2);
  pp.query("X_N2", PeleLM::prob_parm->X_N2);
  pp.query("pertmag_cf", PeleLM::prob_parm->pertmag_cf);
  pp.query("pertmag_jet", PeleLM::prob_parm->pertmag_jet);
  pp.query("jet_purity", PeleLM::prob_parm->jet_purity);
  pp.query("bl_thickness", PeleLM::prob_parm->bl_thickness);
  pp.query("init_time", PeleLM::prob_parm->init_time);
  pp.query("double_jet", PeleLM::prob_parm->double_jet);
  pp.query("jet_dir", PeleLM::prob_parm->jet_dir);
  pp.query("cf_dir", PeleLM::prob_parm->cf_dir);
}
