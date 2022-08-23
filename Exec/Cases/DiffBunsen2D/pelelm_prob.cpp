#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   std::string type;
   
   //Definition of the Geometry
   pp.query("Xmax", prob_parm->Xmax);
   pp.query("Xf", prob_parm->Xf);
   pp.query("Xe", prob_parm->Xe);
   pp.query("Xc", prob_parm->Xc);
   pp.query("Ymax", prob_parm->Ymax);
   pp.query("Zmax", prob_parm->Zmax);
   
   //Definition of the Inflow conditions
   pp.query("Yin", prob_parm->Yin);
   
   pp.query("V_fu", prob_parm->V_fu);
   pp.query("V_ox", prob_parm->V_ox);
   pp.query("V_air", prob_parm->V_air);
   
   pp.query("T_fu", prob_parm->T_fu);
   pp.query("T_ox", prob_parm->T_ox);
   pp.query("T_air", prob_parm->T_air);
   pp.query("T_obst", prob_parm->T_obst);
   
   //Definition of the Initial conditions
   pp.query("P_mean", prob_parm->P_mean);
   
   //Definition of the Ignition zone
   pp.query("do_ignition", prob_parm->do_ignition);
   pp.query("ign_rad", prob_parm->ign_rad);
   pp.query("ign_T", prob_parm->ign_T);
   
   //Definition of the fuel parameters
   pp.query("dilution", prob_parm->dilution);

#ifdef PELE_USE_EFIELD
   pp.query("PhiV_y_hi", PeleLM::prob_parm->phiV_hiy);
   pp.query("PhiV_y_lo", PeleLM::prob_parm->phiV_loy);
#endif
}
