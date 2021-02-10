#include <PeleLM.H>

using namespace amrex;

void PeleLM::WritePlotFile() {
   BL_PROFILE("PeleLM::WritePlotFile()");

   const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep);

   int ncomp = 0;

   

   // Plot MF  
   Vector<MultiFab> mf_plt(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
   }

   Vector<std::string> plt_VarsName;

}
