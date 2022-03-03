#include <PeleLM.H>
#include <pelelm_prob.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   pp.query("P_mean", PeleLM::prob_parm->P_mean);
   pp.query("Zst",    PeleLM::prob_parm->Zst);
   pp.query("T_in",   PeleLM::prob_parm->T_in);
   pp.query("U_b",    PeleLM::prob_parm->U_b);

   auto problo = geom[0].ProbLo();
   auto probhi = geom[0].ProbHi();

   PeleLM::prob_parm->splitx = 0.5 * (problo[0] + probhi[0]);
   PeleLM::prob_parm->midtanh = 0.6 * (problo[0] + probhi[0]);
   PeleLM::prob_parm->widthtanh = 0.05 * (problo[0] + probhi[0]); 

   PeleLM::prob_parm->bathID = N2_ID;  
   PeleLM::prob_parm->fuelID = CH4_ID;  
   PeleLM::prob_parm->oxidID = O2_ID; 
}
