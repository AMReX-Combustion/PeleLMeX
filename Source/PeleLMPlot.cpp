#include <PeleLM.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include "PelePhysics.H"
#include <AMReX_ParmParse.H>

using namespace amrex;

namespace { const std::string level_prefix{"Level_"}; }

void GotoNextLine(std::istream& is)
{
       constexpr std::streamsize bl_ignore_max{100000};
           is.ignore(bl_ignore_max, '\n');
}


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
      // Cons Rate
      ncomp += NUM_SPECIES;
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
      for (int n = 0; n < NUM_SPECIES; n++) {
         plt_VarsName.push_back("I_R("+names[n]+")");
      }
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

      if (m_do_react) {
         MultiFab::Copy(mf_plt[lev], m_leveldatareact[lev]->I_R, 0, cnt, NUM_SPECIES, 0);
         cnt += NUM_SPECIES;

         MultiFab::Copy(mf_plt[lev], m_leveldatareact[lev]->functC, 0, cnt, 1, 0);
         cnt += 1;
      }

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

void PeleLM::WriteHeader(const std::string& name, bool is_checkpoint) const
{
    if(ParallelDescriptor::IOProcessor())
    {
        std::string HeaderFileName(name + "/Header");
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream HeaderFile;

        HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        HeaderFile.open(HeaderFileName.c_str(),
                        std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

        if(!HeaderFile.good()) {
            amrex::FileOpenFailed(HeaderFileName);
        }

        HeaderFile.precision(17);
        if(is_checkpoint) {
            HeaderFile << "Checkpoint version: 1\n";
        } else {
            HeaderFile << "HyperCLaw-V1.1\n";
        }

        HeaderFile << finest_level << "\n";

        // Time stepping controls
        HeaderFile << m_nstep << "\n";
        HeaderFile << m_cur_time << "\n";
        HeaderFile << m_dt << "\n";
        HeaderFile << m_prev_dt << "\n";
        //HeaderFile << m_prev_prev_dt << "\n";

        // Geometry
        for(int i = 0; i < AMREX_SPACEDIM; ++i) {
            HeaderFile << Geom(0).ProbLo(i) << ' ';
        }
        HeaderFile << '\n';

        for(int i = 0; i < AMREX_SPACEDIM; ++i)
            HeaderFile << Geom(0).ProbHi(i) << ' ';
        HeaderFile << '\n';

        // BoxArray
        for(int lev = 0; lev <= finest_level; ++lev)
        {
            boxArray(lev).writeOn(HeaderFile);
            HeaderFile << '\n';
        }
    }
}

void PeleLM::WriteCheckPointFile()
{
   BL_PROFILE("PeleLM::WriteCheckPointFile()");
   
   const std::string& checkpointname = amrex::Concatenate(m_check_file, m_nstep);
   
   if (m_verbose) {
      amrex::Print() << "\n Writting checkpoint file: " << checkpointname << "\n";
   }
   
   amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);
   
   bool is_checkpoint = true;
   WriteHeader(checkpointname, is_checkpoint);
   WriteJobInfo(checkpointname);
   
   for(int lev = 0; lev <= finest_level; ++lev)
   {   
      VisMF::Write(m_leveldata_new[lev]->velocity,
                   amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "velocity"));
   
      VisMF::Write(m_leveldata_new[lev]->gp,
                   amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "gradp"));
   
      VisMF::Write(m_leveldata_new[lev]->press,
                   amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p"));
   
      if (!m_incompressible) {
         VisMF::Write(m_leveldata_new[lev]->density,
                      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "density"));
   
         VisMF::Write(m_leveldata_new[lev]->species,
                      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "species"));
   
         VisMF::Write(m_leveldata_new[lev]->rhoh,
                      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "rhoH"));
   
         VisMF::Write(m_leveldata_new[lev]->temp,
                      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "temp"));
   
         VisMF::Write(m_leveldata_new[lev]->rhoRT,
                      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "rhoRT"));
   
         if (m_has_divu) {
            VisMF::Write(m_leveldata_new[lev]->divu,
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "divU"));
         }
   
         if (m_do_react) {
            VisMF::Write(m_leveldatareact[lev]->I_R,
                         amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "I_R"));
         }
      }    
   }   
}

void PeleLM::ReadCheckPointFile()
{
   BL_PROFILE("PeleLM::ReadCheckPointFile()");

   amrex::Print() << "Restarting from checkpoint " << m_restart_file << "\n";

   Real prob_lo[AMREX_SPACEDIM];
   Real prob_hi[AMREX_SPACEDIM];

   /***************************************************************************
   ** Load header: set up problem domain (including BoxArray)                 *
   **              allocate incflo memory (incflo::AllocateArrays)            *
   **              (by calling MakeNewLevelFromScratch)                       *
   ****************************************************************************/

   std::string File(m_restart_file + "/Header");

   VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

   Vector<char> fileCharPtr;
   ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
   std::string fileCharPtrString(fileCharPtr.dataPtr());
   std::istringstream is(fileCharPtrString, std::istringstream::in);   

   std::string line, word;

   // Start reading from checkpoint file 
   
   // Title line
   std::getline(is, line);

   // Finest level
   is >> finest_level;
   GotoNextLine(is);

   // Step count
   is >> m_nstep;
   GotoNextLine(is);

   // Current time
   is >> m_cur_time;
   GotoNextLine(is);

   // Time step size
   is >> m_dt;
   GotoNextLine(is);

   is >> m_prev_dt;
   GotoNextLine(is);

   //is >> m_prev_prev_dt;
   //GotoNextLine(is);

   // Low coordinates of domain bounding box
   std::getline(is, line);
   {
       std::istringstream lis(line);
       int i = 0;
       while(lis >> word)
       {
           prob_lo[i++] = std::stod(word);
       }
   }

   // High coordinates of domain bounding box
   std::getline(is, line);
   {
       std::istringstream lis(line);
       int i = 0;
       while(lis >> word)
       {
           prob_hi[i++] = std::stod(word);
       }
   }

   // Set up problem domain
   RealBox rb(prob_lo, prob_hi);
   Geometry::ResetDefaultProbDomain(rb);
   for (int lev = 0; lev <= max_level; ++lev) {
       SetGeometry(lev, Geometry(Geom(lev).Domain(), rb, Geom(lev).CoordInt(),
                                 Geom(lev).isPeriodic()));
   }

   for(int lev = 0; lev <= finest_level; ++lev)
   {
       // read in level 'lev' BoxArray from Header
       BoxArray ba;
       ba.readFrom(is);
       GotoNextLine(is);

       // Create distribution mapping
       DistributionMapping dm{ba, ParallelDescriptor::NProcs()};
       MakeNewLevelFromScratch(lev, m_cur_time, ba, dm);
   }

   /***************************************************************************
    * Load fluid data                                                         *
    ***************************************************************************/

   // Load the field data
   for(int lev = 0; lev <= finest_level; ++lev)
   {
      VisMF::Read(m_leveldata_new[lev]->velocity,
             		amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "velocity"));
   
      VisMF::Read(m_leveldata_new[lev]->gp,
             		amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "gradp"));
   
      VisMF::Read(m_leveldata_new[lev]->press,
                   amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "p"));
   
      if (!m_incompressible) {
         VisMF::Read(m_leveldata_new[lev]->density,
                 		amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "density"));
   
         VisMF::Read(m_leveldata_new[lev]->species,
                 		amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "species"));
   
         VisMF::Read(m_leveldata_new[lev]->rhoh,
                 		amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "rhoH"));
   
         VisMF::Read(m_leveldata_new[lev]->temp,
                 		amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "temp"));
   
         VisMF::Read(m_leveldata_new[lev]->rhoRT,
                     amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "rhoRT"));
   
         if (m_has_divu) {
            VisMF::Read(m_leveldata_new[lev]->divu,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "divU"));
         }
   
         if (m_do_react) {
            VisMF::Read(m_leveldatareact[lev]->I_R,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_file, level_prefix, "I_R"));
         }
      }    
   }
	if (m_verbose) {
   	amrex::Print() << "Restart complete" << std::endl;
	}
}

void PeleLM::WriteJobInfo(const std::string& path) const
{
   std::string PrettyLine = std::string(78, '=') + "\n";
   std::string OtherLine = std::string(78, '-') + "\n";
   std::string SkipSpace = std::string(8, ' ');

   if (ParallelDescriptor::IOProcessor()) {
      // job_info file with details about the run
      std::ofstream jobInfoFile;
      std::string FullPathJobInfoFile = path;

      FullPathJobInfoFile += "/PeleLMeX_job_info";
      jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

      // job information
      jobInfoFile << PrettyLine;
      jobInfoFile << " PeleLMeX Job Information\n";
      jobInfoFile << PrettyLine;

      jobInfoFile << "number of MPI processes: " << ParallelDescriptor::NProcs() << "\n";
#ifdef _OPENMP
      jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

      jobInfoFile << "\n\n";

      // build information
      jobInfoFile << PrettyLine;
      jobInfoFile << " Build Information\n";
      jobInfoFile << PrettyLine;

      jobInfoFile << "build date:    " << buildInfoGetBuildDate() << "\n";
      jobInfoFile << "build machine: " << buildInfoGetBuildMachine() << "\n";
      jobInfoFile << "build dir:     " << buildInfoGetBuildDir() << "\n";
      jobInfoFile << "AMReX dir:     " << buildInfoGetAMReXDir() << "\n";

      jobInfoFile << "\n";

      jobInfoFile << "COMP:          " << buildInfoGetComp() << "\n";
      jobInfoFile << "COMP version:  " << buildInfoGetCompVersion() << "\n";
      jobInfoFile << "C++ compiler:  " << buildInfoGetCXXName() << "\n";
      jobInfoFile << "C++ flags:     " << buildInfoGetCXXFlags() << "\n";

      jobInfoFile << "\n";

      const char* githash1 = buildInfoGetGitHash(1);
      const char* githash2 = buildInfoGetGitHash(2);
      const char* githash3 = buildInfoGetGitHash(3);

      if (strlen(githash1) > 0) {
         jobInfoFile << "PeleLMeX     git describe: " << githash1 << "\n";
      }
      if (strlen(githash2) > 0) {
         jobInfoFile << "AMReX        git describe: " << githash2 << "\n";
      }
      if (strlen(githash3) > 0) {
         jobInfoFile << "PelePhysics  git describe: " << githash3 << "\n";
      }

      jobInfoFile << "\n\n";

      // grid information
      jobInfoFile << PrettyLine;
      jobInfoFile << " Grid Information\n";
      jobInfoFile << PrettyLine;

      for(int lev = 0; lev <= finest_level; lev++)
      {
          jobInfoFile << " level: " << lev << "\n";
          jobInfoFile << "   number of boxes = " << grids[lev].size() << "\n";
          jobInfoFile << "   maximum zones   = ";
          for(int idim = 0; idim < AMREX_SPACEDIM; idim++)
          {
              jobInfoFile << geom[lev].Domain().length(idim) << " "; 
          }
          jobInfoFile << "\n\n";
      }

      jobInfoFile << "\n\n";

      // runtime parameters
      jobInfoFile << PrettyLine;
      jobInfoFile << " Inputs File Parameters\n";
      jobInfoFile << PrettyLine;

      ParmParse::dumpTable(jobInfoFile, true);

      jobInfoFile.close();
    }
}


