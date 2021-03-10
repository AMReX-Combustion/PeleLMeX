#include <PeleLM.H>
#include <AMReX_PlotFileUtil.H>
#include <EOS.H>

using namespace amrex;

void PeleLM::WritePlotFile() {
   BL_PROFILE("PeleLM::WritePlotFile()");

   const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep);

   if (m_verbose) {
      amrex::Print() << " Dumping plotfile: " << plotfilename << "\n";
   }

   //----------------------------------------------------------------
   // Number of components
   int ncomp = 0;
   if (m_incompressible) {
      ncomp = 2*AMREX_SPACEDIM + m_derivePlotVarCount;
   } else {
      ncomp = NVAR + AMREX_SPACEDIM + m_derivePlotVarCount;
      if (m_has_divu) {
         ncomp += 1;
      }
   }

   //----------------------------------------------------------------
   // Plot MultiFabs
   Vector<MultiFab> mf_plt(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
   }

   //----------------------------------------------------------------
   // Components names
   Vector<std::string> names;
   EOS::speciesNames(names);

   Vector<std::string> plt_VarsName;
   plt_VarsName.push_back("x_velocity");
#if ( AMREX_SPACEDIM > 1 )
   plt_VarsName.push_back("y_velocity");
#if ( AMREX_SPACEDIM > 2 )
   plt_VarsName.push_back("z_velocity");
#endif
#endif
   if (!m_incompressible) {
      plt_VarsName.push_back("density");
      for (int n = 0; n < NUM_SPECIES; n++) {
         plt_VarsName.push_back("rho.Y("+names[n]+")");
      }
      plt_VarsName.push_back("rhoh");
      plt_VarsName.push_back("temp");
      plt_VarsName.push_back("RhoRT");
      if (m_has_divu) {
         plt_VarsName.push_back("divu");
      }
   }

   plt_VarsName.push_back("gradp_x");
#if ( AMREX_SPACEDIM > 1 )
   plt_VarsName.push_back("gradp_y");
#if ( AMREX_SPACEDIM > 2 )
   plt_VarsName.push_back("gradp_z");
#endif
#endif
   for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
      plt_VarsName.push_back(m_derivePlotVars[ivar]);
   }
   
   //----------------------------------------------------------------
   // Fill the plot MultiFabs
   for (int lev = 0; lev <= finest_level; ++lev) {
      int cnt = 0;
      MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->velocity, 0, cnt, AMREX_SPACEDIM, 0);
      cnt += AMREX_SPACEDIM;
      if (!m_incompressible) {
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->density, 0, cnt, 1, 0);
         cnt += 1;
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->species, 0, cnt, NUM_SPECIES, 0);
         cnt += NUM_SPECIES;
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->rhoh, 0, cnt, 1, 0);
         cnt += 1;
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->temp, 0, cnt, 1, 0);
         cnt += 1;
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->rhoRT, 0, cnt, 1, 0);
         cnt += 1;
         if (m_has_divu) {
            MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->divu, 0, cnt, 1, 0);
            cnt += 1;
         }
      }
      MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->gp, 0, cnt,AMREX_SPACEDIM,0);
      cnt += AMREX_SPACEDIM;

      for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
         std::unique_ptr<MultiFab> mf;
         mf = derive(m_derivePlotVars[ivar], m_cur_time, lev, 0);
         MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, 1, 0);
         cnt += 1;
      }
   }

   // No SubCycling, all levels the same step.
   Vector<int> istep(finest_level + 1, m_nstep);

   amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt),
                                  plt_VarsName, Geom(), m_cur_time, istep, refRatio());

}
