#include <PeleLM.H>
#include <AMReX_ParmParse.H>

void PeleLM::readProbParm()
{
   amrex::ParmParse pp("prob");
   
   std::string type;
   pp.query("P_mean",   PeleLM::prob_parm->P_mean);
   pp.query("init_T", PeleLM::prob_parm->T0);
   pp.query("init_vel", PeleLM::prob_parm->vel);
   pp.query("init_N2", PeleLM::prob_parm->Y_N2);
   pp.query("init_O2", PeleLM::prob_parm->Y_O2);
   // Find the number of redistributions during particle initialization
   pp.query("init_redist", PeleLM::prob_parm->numRedist);
   pp.query("num_particles", PeleLM::prob_parm->partNum);
   std::array<amrex::Real, AMREX_SPACEDIM> pvel;
   pp.query<amrex::Real>("part_vel", pvel);
   for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
     PeleLM::prob_parm->partVel[dir] = pvel[dir];
   }
   pp.get("part_dia", PeleLM::prob_parm->partDia);
   pp.get("part_temp", PeleLM::prob_parm->partTemp);
}
