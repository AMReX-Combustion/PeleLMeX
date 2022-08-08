#include <PeleLM.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include "PelePhysics.H"
#include <PltFileManager.H>
#include <AMReX_ParmParse.H>
#include <PeleLMBCfill.H>
#include <AMReX_FillPatchUtil.H>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

#ifdef PELELM_USE_SPRAY
#include "SprayParticles.H"
#endif
#ifdef PELELM_USE_SOOT
#include "SootModel.H"
#endif
using namespace amrex;

namespace { const std::string level_prefix{"Level_"}; }

void GotoNextLine(std::istream& is)
{
       constexpr std::streamsize bl_ignore_max{100000};
           is.ignore(bl_ignore_max, '\n');
}

void PeleLM::WriteDebugPlotFile(const Vector<const MultiFab*> &a_MF,
                                const std::string &pltname)
{
   int nComp = a_MF[0]->nComp();
   Vector<std::string> names(nComp);
   for (int n = 0; n < nComp; n++) {
      names[n] = "comp"+std::to_string(n);
   }
   Vector<int> istep(finest_level + 1, m_nstep);
#ifdef AMREX_USE_HDF5
   if (m_write_hdf5_pltfile) {
       amrex::WriteMultiLevelPlotfileHDF5(pltname, finest_level + 1, a_MF,
                                          names, Geom(), m_cur_time, istep, refRatio());
   } else
#endif
   {
       amrex::WriteMultiLevelPlotfile(pltname, finest_level + 1, a_MF,
                                      names, Geom(), m_cur_time, istep, refRatio());
   }
}

void PeleLM::WritePlotFile() {
   BL_PROFILE("PeleLM::WritePlotFile()");

   const std::string& plotfilename = amrex::Concatenate(m_plot_file, m_nstep, m_ioDigits);

   if (m_verbose) {
      amrex::Print() << " Dumping plotfile: " << plotfilename << "\n";
   }

   //----------------------------------------------------------------
   // Number of components
   int ncomp = 0;

   // State
   if (m_incompressible) {
      // Velocity + pressure gradients
      ncomp = 2*AMREX_SPACEDIM;
   } else {
      // State + pressure gradients
      if (m_plot_grad_p) {
         ncomp = NVAR + AMREX_SPACEDIM;
      } else {
         ncomp = NVAR;
      }
      // Make the plot lighter by dropping species by default
      if (!m_plotStateSpec) ncomp -= NUM_SPECIES;
      if (m_has_divu) {
         ncomp += 1;
      }
   }

   // Reactions
   if (m_do_react && !m_skipInstantRR && m_plot_react) {
      // Cons Rate
      ncomp += nCompIR();
      // FunctCall
      ncomp += 1;
      // Extras:
      if (m_plotHeatRelease) ncomp += 1;
      //if (m_plotChemDiag) ncomp += 1;     // TODO
   }

#ifdef AMREX_USE_EB
   // Include volume fraction in plotfile
   ncomp += 1;
#endif

   // Derive
   int deriveEntryCount = 0;
   for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
      const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
      deriveEntryCount += rec->numDerive();
   }
   ncomp += deriveEntryCount;
#ifdef PELELM_USE_SPRAY
   if (do_spray_particles) {
     ncomp += spray_derive_vars.size();
   }
#endif

#ifdef PELE_USE_EFIELD
   if (m_do_extraEFdiags) {
       ncomp += NUM_IONS * AMREX_SPACEDIM;
   }
#endif

   //----------------------------------------------------------------
   // Plot MultiFabs
   Vector<MultiFab> mf_plt(finest_level + 1);
   for (int lev = 0; lev <= finest_level; ++lev) {
      mf_plt[lev].define(grids[lev], dmap[lev], ncomp, 0, MFInfo(), Factory(lev));
   }

   //----------------------------------------------------------------
   // Components names
   Vector<std::string> names;
   pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(names);

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
      if (m_plotStateSpec) {
         for (int n = 0; n < NUM_SPECIES; n++) {
            plt_VarsName.push_back("rho.Y("+names[n]+")");
         }
      }
      plt_VarsName.push_back("rhoh");
      plt_VarsName.push_back("temp");
      plt_VarsName.push_back("RhoRT");
#ifdef PELE_USE_EFIELD
      plt_VarsName.push_back("nE");
      plt_VarsName.push_back("phiV");
#endif
#ifdef PELELM_USE_SOOT
      for (int mom = 0; mom < NUMSOOTVAR; mom++) {
        std::string sootname = soot_model->sootVariableName(mom);
        plt_VarsName.push_back(sootname);
      }
#endif
      if (m_has_divu) {
         plt_VarsName.push_back("divu");
      }
   }

   if (m_plot_grad_p) {
      plt_VarsName.push_back("gradpx");
#if ( AMREX_SPACEDIM > 1 )
      plt_VarsName.push_back("gradpy");
#if ( AMREX_SPACEDIM > 2 )
      plt_VarsName.push_back("gradpz");
#endif
#endif
   }

   if (m_do_react  && !m_skipInstantRR && m_plot_react) {
      for (int n = 0; n < NUM_SPECIES; n++) {
         plt_VarsName.push_back("I_R("+names[n]+")");
      }
#ifdef PELE_USE_EFIELD
      plt_VarsName.push_back("I_R(nE)");
#endif
      plt_VarsName.push_back("FunctCall");
      // Extras:
      if (m_plotHeatRelease) plt_VarsName.push_back("HeatRelease");
      //if (m_plotChemDiag) ncomp += 1;     // TODO
   }

#ifdef AMREX_USE_EB
   plt_VarsName.push_back("volFrac");
#endif

   for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
      const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
      for (int dvar = 0; dvar < rec->numDerive(); dvar++ ) {
         plt_VarsName.push_back(rec->variableName(dvar));
      }
   }
#ifdef PELELM_USE_SPRAY
   if (spray_derive_vars.size() > 0) {
     // We need virtual particles for the lower levels
     setupVirtualParticles(0);
     for (int ivar = 0; ivar < spray_derive_vars.size(); ivar++) {
       plt_VarsName.push_back(spray_derive_vars[ivar]);
     }
   }
#endif

#ifdef PELE_USE_EFIELD
   if (m_do_extraEFdiags) {
       for (int ivar = 0; ivar < NUM_IONS; ++ivar) {
           for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
               std::string dir = (idim == 0) ? "X" : ( (idim == 1) ? "Y" : "Z");
               plt_VarsName.push_back("DriftFlux_"+names[NUM_SPECIES-NUM_IONS+ivar]+"_"+dir);
           }
       }
   }
#endif

   //----------------------------------------------------------------
   // Fill the plot MultiFabs
   for (int lev = 0; lev <= finest_level; ++lev) {
      int cnt = 0;
      if (m_incompressible) {
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, 0, cnt, AMREX_SPACEDIM, 0);
         cnt += AMREX_SPACEDIM;
      } else {
         // Velocity and density
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, 0, cnt, AMREX_SPACEDIM+1, 0);
         cnt += AMREX_SPACEDIM+1;
         // Species only if requested
         if (m_plotStateSpec) {
            MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, FIRSTSPEC, cnt, NUM_SPECIES, 0);
            cnt += NUM_SPECIES;
         }
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, RHOH, cnt, 3, 0);
         cnt += 3;
#ifdef PELE_USE_EFIELD
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, NE, cnt, 2, 0);
         cnt += 2;
#endif
#ifdef PELELM_USE_SOOT
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->state, FIRSTSOOT, cnt, NUMSOOTVAR, 0);
         cnt += NUMSOOTVAR;
#endif
         if (m_has_divu) {
            MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->divu, 0, cnt, 1, 0);
            cnt += 1;
         }
      }
      if (m_plot_grad_p) {
         MultiFab::Copy(mf_plt[lev], m_leveldata_new[lev]->gp, 0, cnt,AMREX_SPACEDIM,0);
         cnt += AMREX_SPACEDIM;
      }

      if (m_do_react  && !m_skipInstantRR && m_plot_react) {
         MultiFab::Copy(mf_plt[lev], m_leveldatareact[lev]->I_R, 0, cnt, nCompIR(), 0);
         cnt += nCompIR();

         MultiFab::Copy(mf_plt[lev], m_leveldatareact[lev]->functC, 0, cnt, 1, 0);
         cnt += 1;

         if (m_plotHeatRelease) {
            std::unique_ptr<MultiFab> mf;
            mf.reset( new MultiFab(grids[lev],dmap[lev],1,0));
            getHeatRelease(lev, mf.get());
            MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, 1, 0);
            cnt += 1;
         }
      }

#ifdef AMREX_USE_EB
      MultiFab::Copy(mf_plt[lev], EBFactory(lev).getVolFrac(), 0, cnt, 1, 0);
      cnt += 1;
#endif

      for (int ivar = 0; ivar < m_derivePlotVarCount; ivar++ ) {
         std::unique_ptr<MultiFab> mf;
         mf = derive(m_derivePlotVars[ivar], m_cur_time, lev, 0);
         MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, mf->nComp(), 0);
         cnt += mf->nComp();
      }
#ifdef PELELM_USE_SPRAY
      if (spray_derive_vars.size() > 0) {
        int num_spray_derive = spray_derive_vars.size();
        mf_plt[lev].setVal(0., cnt, num_spray_derive);
        theSprayPC()->computeDerivedVars(
          mf_plt[lev], lev, cnt, spray_derive_vars, spray_fuel_names);
        if (lev < finest_level) {
          MultiFab tmp_plt(grids[lev], dmap[lev], num_spray_derive, 0, MFInfo(), Factory(lev));
          tmp_plt.setVal(0.);
          theVirtPC()->computeDerivedVars(
            tmp_plt, lev, 0, spray_derive_vars, spray_fuel_names);
          MultiFab::Add(mf_plt[lev], tmp_plt, 0, cnt, num_spray_derive, 0);
        }
        cnt += num_spray_derive;
      }
#endif
#ifdef PELE_USE_EFIELD
      if (m_do_extraEFdiags) {
          MultiFab::Copy(mf_plt[lev], *m_ionsFluxes[lev], 0, cnt, m_ionsFluxes[lev]->nComp(),0);
      }
#endif
#ifdef AMREX_USE_EB
      EB_set_covered(mf_plt[lev],0.0);
#endif
   }


   // No SubCycling, all levels the same step.
   Vector<int> istep(finest_level + 1, m_nstep);

#ifdef AMREX_USE_HDF5
   if (m_write_hdf5_pltfile) {
       amrex::WriteMultiLevelPlotfileHDF5(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt),
                                          plt_VarsName, Geom(), m_cur_time, istep, refRatio());
   } else
#endif
   {
       amrex::WriteMultiLevelPlotfile(plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt),
                                      plt_VarsName, Geom(), m_cur_time, istep, refRatio());
   }

#ifdef PELELM_USE_SPRAY
   if (theSprayPC() != nullptr && do_spray_particles) {
     bool is_spraycheck = false;
     for (int lev = 0; lev <= finest_level; ++lev) {
       theSprayPC()->SprayParticleIO(
         lev, is_spraycheck, write_spray_ascii_files, plotfilename, spray_fuel_names);
       // Remove virtual particles that were made for derived variables
       removeVirtualParticles(lev);
     }
   }
#endif
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
        
        HeaderFile << m_nstep << "\n";

#ifdef AMREX_USE_EB
        HeaderFile << m_EB_generate_max_level << "\n";
#endif

        HeaderFile << m_cur_time << "\n";
        HeaderFile << m_dt << "\n";
        HeaderFile << m_prev_dt << "\n";

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

        // Ambient pressure and typvals
        HeaderFile << m_pNew << "\n";
        for (int n = 0; n < typical_values.size(); n++) {
            HeaderFile << typical_values[n] << "\n";
        }
    }
}

void PeleLM::WriteCheckPointFile()
{
   BL_PROFILE("PeleLM::WriteCheckPointFile()");

   const std::string& checkpointname = amrex::Concatenate(m_check_file, m_nstep, m_ioDigits);

   if (m_verbose) {
      amrex::Print() << "\n Writting checkpoint file: " << checkpointname << "\n";
   }

   amrex::PreBuildDirectorHierarchy(checkpointname, level_prefix, finest_level + 1, true);

   bool is_checkpoint = true;
   WriteHeader(checkpointname, is_checkpoint);
   WriteJobInfo(checkpointname);

   for(int lev = 0; lev <= finest_level; ++lev)
   {
      VisMF::Write(m_leveldata_new[lev]->state,
                   amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "state"));

      VisMF::Write(m_leveldata_new[lev]->gp,
                   amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "gradp"));

      VisMF::Write(m_leveldata_new[lev]->press,
                   amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p"));

      if (!m_incompressible) {
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
#ifdef PELELM_USE_SPRAY
   if (theSprayPC() != nullptr && do_spray_particles) {
     int write_ascii = 0; // Not for checkpoints
     bool is_spraycheck = true;
     for (int lev = 0; lev <= finest_level; ++lev) {
       theSprayPC()->SprayParticleIO(
         lev, is_spraycheck, write_ascii, checkpointname, spray_fuel_names);
     }
   }
#endif
}

void PeleLM::ReadCheckPointFile()
{
   BL_PROFILE("PeleLM::ReadCheckPointFile()");

   amrex::Print() << "Restarting from checkpoint " << m_restart_chkfile << "\n";

   Real prob_lo[AMREX_SPACEDIM];
   Real prob_hi[AMREX_SPACEDIM];

   /***************************************************************************
   ** Load header: set up problem domain (including BoxArray)                 *
   **              allocate PeleLM memory (PeleLM::AllocateArrays)            *
   **              (by calling MakeNewLevelFromScratch)                       *
   ****************************************************************************/

   std::string File(m_restart_chkfile + "/Header");

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
   int chk_finest_level = 0;
   is >> chk_finest_level;
   GotoNextLine(is);
   finest_level = std::min(chk_finest_level,max_level);

   // Step count
   is >> m_nstep;
   GotoNextLine(is);

#ifdef AMREX_USE_EB
   // Finest level at which EB was generated
   // actually used independently, so just skip ...
   std::getline(is, line);

   // ... but to be backward compatible, if we get a float,
   // let's assume it's m_cur_time
   if (line.find('.') != std::string::npos) {
      m_cur_time = std::stod(line);
   } else {
      // Skip line and read current time
      is >> m_cur_time;
      GotoNextLine(is);
   }
#else 

   // Current time
   is >> m_cur_time;
   GotoNextLine(is);
#endif

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

   // deal with typval and P_amb
   is >> m_pNew;
   GotoNextLine(is);
   m_pOld = m_pNew;
   for (int n = 0; n < typical_values.size(); n++) {
       is >> typical_values[n];
       GotoNextLine(is);
   }

   /***************************************************************************
    * Load fluid data                                                         *
    ***************************************************************************/

   // Load the field data
   for(int lev = 0; lev <= finest_level; ++lev)
   {
#ifdef PELE_USE_EFIELD
      if (!m_restart_nonEF) {
         VisMF::Read(m_leveldata_new[lev]->state,
                     amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "state"));
      } else {
         // The chk state is 2 component shorter since phiV and nE aren't in it
         MultiFab stateTemp(grids[lev],dmap[lev],NVAR-2,m_nGrowState);
         VisMF::Read(stateTemp,
                     amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "state"));
         MultiFab::Copy(m_leveldata_new[lev]->state, stateTemp, 0, 0, NVAR-2, m_nGrowState);
      }
#else
      VisMF::Read(m_leveldata_new[lev]->state,
                  amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "state"));
#endif

      VisMF::Read(m_leveldata_new[lev]->gp,
                  amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "gradp"));

      VisMF::Read(m_leveldata_new[lev]->press,
                  amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "p"));

      if (!m_incompressible) {
         if (m_has_divu) {
            VisMF::Read(m_leveldata_new[lev]->divu,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "divU"));
         }

#ifdef PELE_USE_EFIELD
         if (!m_restart_nonEF) {
            if (m_do_react) {
               VisMF::Read(m_leveldatareact[lev]->I_R,
                           amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "I_R"));
            }
         } else {
            // I_R for non-EF simulation is one component shorted, need to account for that.
            if (m_do_react) {
               MultiFab I_Rtemp(grids[lev],dmap[lev],NUM_SPECIES,0);
               VisMF::Read(I_Rtemp,amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "I_R"));
               MultiFab::Copy(m_leveldatareact[lev]->I_R,I_Rtemp,0,0,NUM_SPECIES,0);
            }

            // Initialize nE & phiV
            m_leveldata_new[lev]->state.setVal(0.0,NE,2,m_nGrowState);
         }
#else
         if (m_do_react) {
            VisMF::Read(m_leveldatareact[lev]->I_R,
                        amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "I_R"));
         }
#endif
      }
   }
   if (m_verbose) {
      amrex::Print() << "Restart complete" << std::endl;
   }
}

void PeleLM::initLevelDataFromPlt(int a_lev,
                                  const std::string &a_dataPltFile)
{
   if (m_incompressible) {
      Abort(" initializing data from a pltfile only available for low-Mach simulations");
   }

   amrex::Print() << " initData on level " << a_lev << " from pltfile " << a_dataPltFile << "\n";
   if(pltfileSource == "LM"){
     amrex::Print() << " Assuming pltfile was generated in LM/LMeX \n"; 
   }
   else if(pltfileSource == "C"){
     amrex::Print() << " Assuming pltfile was generated in PeleC \n"; 
   }

   // Use PelePhysics PltFileManager
   pele::physics::pltfilemanager::PltFileManager pltData(a_dataPltFile);
   Vector<std::string> plt_vars = pltData.getVariableList();

   // Find required data in pltfile
   Vector<std::string> spec_names;
   pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(spec_names);
   int idT = -1, idV = -1, idY = -1, nSpecPlt = 0;
#ifdef PELE_USE_EFIELD
   int inE = -1, iPhiV = -1; 
#endif
   for (int i = 0; i < plt_vars.size(); ++i) {
      std::string firstChars = plt_vars[i].substr(0, 2);
      
      if(pltfileSource == "LM"){
        if (plt_vars[i] == "temp")            idT = i; 
      }
      else if(pltfileSource == "C"){
        if (plt_vars[i] == "Temp")            idT = i; 
      }
      if (plt_vars[i] == "x_velocity")      idV = i; 
      if (plt_vars[i] == "x_velocity")      idV = i;
      if (firstChars == "Y(" && idY < 0 ) {  // species might not be ordered in the order of the current mech.
         idY = i;
      }
      if (firstChars == "Y(")                nSpecPlt += 1;
#ifdef PELE_USE_EFIELD
      if (plt_vars[i] == "nE")              inE = i;
      if (plt_vars[i] == "phiV")            iPhiV = i;
#endif
   }
   if ( idY < 0 ) {
      Abort("Coudn't find species mass fractions in pltfile");
   }
   else if(idT < 0){
      Abort("Coudn't find temperature in pltfile");
   }
   Print() << " " << nSpecPlt << " species found in pltfile, starting with " << plt_vars[idY] << "\n";

   // Get level data
   auto ldata_p = getLevelDataPtr(a_lev,AmrNewTime);


   // Velocity
   pltData.fillPatchFromPlt(a_lev, geom[a_lev], idV, VELX, AMREX_SPACEDIM,
                            ldata_p->state);

   // Temperature
   pltData.fillPatchFromPlt(a_lev, geom[a_lev], idT, TEMP, 1,
                            ldata_p->state);

   // Species
   // Hold the species in temporary MF before copying to level data
   // in case the number of species differs.
   MultiFab speciesPlt(grids[a_lev], dmap[a_lev], nSpecPlt, 0);
   pltData.fillPatchFromPlt(a_lev, geom[a_lev], idY, 0, nSpecPlt,
                            speciesPlt);
   for (int i = 0; i < NUM_SPECIES; i++) {
      std::string specString = "Y("+spec_names[i]+")";
      int foundSpec = 0;
      for (int iplt = 0; iplt < nSpecPlt; iplt++) {
         if ( specString == plt_vars[idY+iplt] ) {
            MultiFab::Copy(ldata_p->state, speciesPlt, iplt, FIRSTSPEC+i, 1, 0);
            foundSpec = 1;
         }
      }
      if (!foundSpec) ldata_p->state.setVal(0.0,FIRSTSPEC+i,1);
   }

   // Converting units when pltfile is coming from PeleC solution
   if(pltfileSource == "C"){
      amrex::Print() << " Converting CGS to MKS units... \n";
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
         const Box& bx = mfi.tilebox();
         auto  const &vel_arr   = ldata_p->state.array(mfi,VELX);
         amrex::ParallelFor(bx, [=]
         AMREX_GPU_DEVICE (int i, int j, int k) noexcept
         {
            for (int n = 0; n < AMREX_SPACEDIM; n++){
               amrex::Real vel_mks = vel_arr(i,j,k,n) * 0.01;
               vel_arr(i,j,k,n) = vel_mks;
            }
         });
      }
   }

#ifdef PELE_USE_EFIELD
   // nE
   pltData.fillPatchFromPlt(a_lev, geom[a_lev], inE, NE, 1,
                            ldata_p->state);
   // phiV
   pltData.fillPatchFromPlt(a_lev, geom[a_lev], iPhiV, PHIV, 1,
                            ldata_p->state);
#endif

   // Pressure and pressure gradients to zero
   ldata_p->press.setVal(0.0);
   ldata_p->gp.setVal(0.0);

   ProbParm const* lprobparm = prob_parm_d;

   // Enforce rho and rhoH consistent with temperature and mixture
   // TODO the above handles species mapping (to some extent), but nothing enforce
   // sum of Ys = 1
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
   for (MFIter mfi(ldata_p->state,TilingIfNotGPU()); mfi.isValid(); ++mfi)
   {
      const Box& bx = mfi.tilebox();
      auto  const &rho_arr   = ldata_p->state.array(mfi,DENSITY);
      auto  const &rhoY_arr  = ldata_p->state.array(mfi,FIRSTSPEC);
      auto  const &rhoH_arr  = ldata_p->state.array(mfi,RHOH);
      auto  const &temp_arr  = ldata_p->state.array(mfi,TEMP);
      amrex::ParallelFor(bx, [=]
      AMREX_GPU_DEVICE (int i, int j, int k) noexcept
      {
          auto eos = pele::physics::PhysicsType::eos();
          Real massfrac[NUM_SPECIES] = {0.0};
          Real sumYs = 0.0;
          for (int n = 0; n < NUM_SPECIES; n++){
             massfrac[n] = rhoY_arr(i,j,k,n);
             if (n != N2_ID) {
                sumYs += massfrac[n];
             }
          }
          massfrac[N2_ID] = 1.0 - sumYs;

          // Get density
          Real P_cgs = lprobparm->P_mean * 10.0;
          Real rho_cgs = 0.0;
          eos.PYT2R(P_cgs, massfrac, temp_arr(i,j,k), rho_cgs);
          rho_arr(i,j,k) = rho_cgs * 1.0e3;

          // Get enthalpy
          Real h_cgs = 0.0;
          eos.TY2H(temp_arr(i,j,k), massfrac, h_cgs);
          rhoH_arr(i,j,k) = h_cgs * 1.0e-4 * rho_arr(i,j,k);

          // Fill rhoYs
          for (int n = 0; n < NUM_SPECIES; n++){
             rhoY_arr(i,j,k,n) = massfrac[n] * rho_arr(i,j,k);
          }
      });
   }

   // Initialize thermodynamic pressure
   setThermoPress(a_lev, AmrNewTime);
   if (m_has_divu) {
      ldata_p->divu.setVal(0.0);
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
#ifdef AMREX_USE_OMP
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
      const char* githash4 = buildInfoGetGitHash(4);

      if (strlen(githash1) > 0) {
         jobInfoFile << "PeleLMeX     git describe: " << githash1 << "\n";
      }
      if (strlen(githash2) > 0) {
         jobInfoFile << "AMReX        git describe: " << githash2 << "\n";
      }
      if (strlen(githash3) > 0) {
         jobInfoFile << "PelePhysics  git describe: " << githash3 << "\n";
      }
      if (strlen(githash4) > 0) {
         jobInfoFile << "AMREX-Hydro  git describe: " << githash3 << "\n";
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
