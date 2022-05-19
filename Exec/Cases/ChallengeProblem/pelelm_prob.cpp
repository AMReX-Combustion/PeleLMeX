#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   // Chamber conditions
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("T_mean", PeleLM::prob_parm->T_mean);
   pp.query("Y_CH4_chamber", PeleLM::prob_parm->Y_CH4_chamber);
   pp.query("Y_O2_chamber", PeleLM::prob_parm->Y_O2_chamber);

   // Injection parameters
   pp.query("nholes", PeleLM::prob_parm->nholes);
   pp.query("cone_angle", PeleLM::prob_parm->cone_angle);
   pp.query("centx", PeleLM::prob_parm->centx);
   pp.query("centy", PeleLM::prob_parm->centy);
   pp.query("r_circ", PeleLM::prob_parm->r_circ);
   pp.query("r_hole", PeleLM::prob_parm->r_hole);
   pp.query("T_jet", PeleLM::prob_parm->T_jet);
   pp.query("vel_jet", PeleLM::prob_parm->vel_jet);
   pp.query("injection_start", PeleLM::prob_parm->inj_start);
   pp.query("injection_duration", PeleLM::prob_parm->inj_dur);
   pp.query("tau", PeleLM::prob_parm->tau);
   pp.query("Z", PeleLM::prob_parm->Z);
}
