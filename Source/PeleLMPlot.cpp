#include <PeleLM.H>
#include <AMReX_PlotFileUtil.H>
#include "PelePhysics.H"

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

   // State
   if (m_incompressible) {
      ncomp = 2*AMREX_SPACEDIM;
   } else {
      ncomp = NVAR + AMREX_SPACEDIM;
      if (m_has_divu) {
         ncomp += 1;
      }
   }

   // Reactions
   if (m_do_react) {
      // FunctCall
      ncomp += 1;
   }

   // Derive
   int deriveEntryCount = 0;
   for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
      const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
      deriveEntryCount += rec->numDerive();
   }
   ncomp += deriveEntryCount;

   //----------------------------------------------------------------
   // Plot MultiFabs
   Vector<MultiFab> mf_plt(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
   }

   //----------------------------------------------------------------
   // Components names
   Vector<std::string> names;
   pele::physics::eos::speciesNames(names);

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

   if (m_do_react) {
      plt_VarsName.push_back("FunctCall");
   }

   for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
      const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
      for (int dvar = 0; dvar < rec->numDerive(); dvar++ ) {
         plt_VarsName.push_back(rec->variableName(dvar));
      }
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

      MultiFab::Copy(mf_plt[lev], m_leveldatareact[lev]->functC, 0, cnt, 1, 0);
      cnt += 1;

      for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
         std::unique_ptr<MultiFab> mf;
         mf = derive(m_derivePlotVars[ivar], m_cur_time, lev, 0);
         MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, mf->nComp(), 0);
         cnt += mf->nComp();
      }
   }

   // No SubCycling, all levels the same step.
   Vector<int> istep(finest_level + 1, m_nstep);

   amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt),
                                  plt_VarsName, Geom(), m_cur_time, istep, refRatio());

}
