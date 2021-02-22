#include <PeleLM.H>
#include <AMReX_PlotFileUtil.H>
#include <EOS.H>

using namespace amrex;

void PeleLM::WritePlotFile() {
   BL_PROFILE("PeleLM::WritePlotFile()");

   const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep);

   int ncomp = 0;
   if (m_incompressible) {
      ncomp = AMREX_SPACEDIM;
   } else {
      ncomp = NVAR;
   }

   // Plot MF  
   Vector<MultiFab> mf_plt(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
   }

   Vector<std::string> names;
   EOS::speciesNames(names);

   Vector<std::string> plt_VarsName;
   plt_VarsName.push_back("U");
#if ( AMREX_SPACEDIM > 1 )
   plt_VarsName.push_back("V");
#if ( AMREX_SPACEDIM > 2 )
   plt_VarsName.push_back("W");
#endif
#endif
   if (!m_incompressible) {
      plt_VarsName.push_back("density");
      for (int n = 0; n < NUM_SPECIES; n++) {
         plt_VarsName.push_back("rho.Y("+names[n]+")");
      }
      plt_VarsName.push_back("rhoH");
      plt_VarsName.push_back("temp");
      plt_VarsName.push_back("rhoRT");
   }
   
   for (int lev = 0; lev <= finest_level; ++lev) {
      MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->velocity, 0, VELX, AMREX_SPACEDIM, 0);
      if (!m_incompressible) {
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->density, 0, DENSITY, 1, 0);
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->species, 0, FIRSTSPEC, NUM_SPECIES, 0);
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->rhoh, 0, RHOH, 1, 0);
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->temp, 0, TEMP, 1, 0);
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->rhoRT, 0, RHORT, 1, 0);
      }
   }

   // No SubCycling, all levels the same step.
   Vector<int> istep(finest_level + 1, m_nstep);

   amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt),
                                  plt_VarsName, Geom(), m_cur_time, istep, refRatio());

}
