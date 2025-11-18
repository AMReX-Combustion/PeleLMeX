#include <PeleLMeX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_buildInfo.H>
#include "PelePhysics.H"
#include <PltFileManager.H>
#include <AMReX_ParmParse.H>
#include <PeleLMeX_BCfill.H>
#include <AMReX_FillPatchUtil.H>
#include <memory>
#ifdef AMREX_USE_EB
#include <AMReX_EBInterpolater.H>
#endif

#ifdef PELE_USE_SPRAY
#include "SprayParticles.H"
#endif
#ifdef PELE_USE_SOOT
#include "SootModel.H"
#endif
#ifdef PELE_USE_RADIATION
#include "PeleLMRad.H"
#endif

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

namespace {
const std::string level_prefix{"Level_"};
}

void
GotoNextLine(std::istream& is)
{
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}

void
PeleLM::WriteDebugPlotFile(
  const amrex::Vector<const amrex::MultiFab*>& a_MF, const std::string& pltname)
{
  const int nComp = a_MF[0]->nComp();
  amrex::Vector<std::string> names(nComp);
  for (int n = 0; n < nComp; ++n) {
    names[n] = "comp" + std::to_string(n);
  }
  amrex::Vector<int> istep(finest_level + 1, m_nstep);
#ifdef AMREX_USE_HDF5
  if (m_write_hdf5_pltfile) {
    amrex::WriteMultiLevelPlotfileHDF5(
      pltname, finest_level + 1, a_MF, names, Geom(), m_cur_time, istep,
      refRatio());
  } else
#endif
  {
    amrex::WriteMultiLevelPlotfile(
      pltname, finest_level + 1, a_MF, names, Geom(), m_cur_time, istep,
      refRatio());
  }
}

void
PeleLM::WritePlotFile()
{
  BL_PROFILE("PeleLMeX::WritePlotFile()");

  const std::string& plotfilename =
    amrex::Concatenate(m_plot_file, m_nstep, m_ioDigits);

  if (m_verbose != 0) {
    amrex::Print() << "\n Writing plotfile: " << plotfilename << "\n";
  }

  //----------------------------------------------------------------
  // Delete plotfiles if present and requested (and have same name)
  if (m_plot_overwrite) {
    if (amrex::ParallelContext::IOProcessorSub()) {
      if (amrex::FileExists(plotfilename)) {
        amrex::FileSystem::RemoveAll(plotfilename);
      }
    }
  }

  amrex::VisMF::SetNOutFiles(m_nfiles);

  //----------------------------------------------------------------
  // Average down the state
  averageDownState(AmrNewTime);
  if (m_nAux > 0) {
    averageDownAux(AmrNewTime);
  }
  // Get consistent reaction data across level
  if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
    averageDownReaction();
  }

  //----------------------------------------------------------------
  // Number of components
  int ncomp = 0;

  // State
  if (m_incompressible != 0) {
    // Velocity + pressure gradients
    ncomp = 2 * AMREX_SPACEDIM;
  } else {
    // State + pressure gradients
    if (m_plot_grad_p != 0) {
      ncomp = NVAR + AMREX_SPACEDIM;
    } else {
      ncomp = NVAR;
    }
    // Make the plot lighter by dropping species by default
    if (m_plotStateSpec == 0) {
      ncomp -= NUM_SPECIES;
    }
    if (m_has_divu != 0) {
      ncomp += 1;
    }
  }

  ncomp += m_nAux;

  // Reactions
  if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
    // Cons Rate
    ncomp += nCompIR();
    // FunctCall
    ncomp += 1;
    // Extras:
    if (m_plotHeatRelease != 0) {
      ncomp += 1;
    }
  }

#ifdef AMREX_USE_EB
  // Include volume fraction in plotfile
  ncomp += 1;
#endif

#ifdef PELE_USE_RADIATION
  if (do_rad_solve) {
    ncomp += 3;
  }
#endif

  // Derive
  int deriveEntryCount = 0;
  for (int ivar = 0; ivar < m_derivePlotVarCount; ++ivar) {
    const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
    deriveEntryCount += rec->numDerive();
  }
  ncomp += deriveEntryCount;
#ifdef PELE_USE_SPRAY
  if (do_spray_particles) {
    ncomp += SprayParticleContainer::NumDeriveVars();
    if (SprayParticleContainer::plot_spray_src) {
      ncomp += AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;
    }
  }
#endif

#ifdef PELE_USE_PLASMA
  if (m_do_extraEFdiags) {
    ncomp += NUM_IONS * AMREX_SPACEDIM;
  }
#endif

  if (m_do_les && m_plot_les) {
    ncomp += 1;
  }

  if (m_plot_extSource) {
    ncomp += NVAR;
  }

  //----------------------------------------------------------------
  // Plot MultiFabs
  amrex::Vector<amrex::MultiFab> mf_plt;
  mf_plt.reserve(finest_level + 1);
  for (int lev = 0; lev <= finest_level; ++lev) {
    mf_plt.emplace_back(
      grids[lev], dmap[lev], ncomp, 0, amrex::MFInfo(), Factory(lev));
  }

  //----------------------------------------------------------------
  // Components names
  amrex::Vector<std::string> names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    names, &(eos_parms.host_parm()));

  amrex::Vector<std::string> plt_VarsName;
  AMREX_D_TERM(plt_VarsName.push_back("x_velocity");
               , plt_VarsName.push_back("y_velocity");
               , plt_VarsName.push_back("z_velocity"));
  if (m_incompressible == 0) {
    plt_VarsName.push_back("density");
    if (m_plotStateSpec != 0) {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        plt_VarsName.push_back("rho.Y(" + names[n] + ")");
      }
    }
    plt_VarsName.push_back("rhoh");
    plt_VarsName.push_back("temp");
    plt_VarsName.push_back("RhoRT");
#ifdef PELE_USE_PLASMA
    plt_VarsName.push_back("nE");
    plt_VarsName.push_back("phiV");
#endif
#ifdef PELE_USE_SOOT
    for (int mom = 0; mom < NUMSOOTVAR; ++mom) {
      const std::string sootname = soot_model->sootVariableName(mom);
      plt_VarsName.push_back(sootname);
    }
#endif
#ifdef PELE_USE_RADIATION
    if (do_rad_solve) {
      plt_VarsName.push_back("rad.G");
      plt_VarsName.push_back("rad.kappa");
      plt_VarsName.push_back("rad.emis");
    }
#endif
    if (m_has_divu != 0) {
      plt_VarsName.push_back("divu");
    }
  }

  if (m_plot_grad_p != 0) {
    AMREX_D_TERM(plt_VarsName.push_back("gradpx");
                 , plt_VarsName.push_back("gradpy");
                 , plt_VarsName.push_back("gradpz"));
  }

  for (int n = 0; n < m_nAux; ++n) {
    plt_VarsName.push_back(m_aux_names[n]);
  }

  if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
    for (int n = 0; n < NUM_SPECIES; ++n) {
      plt_VarsName.push_back("I_R(" + names[n] + ")");
    }
#ifdef PELE_USE_PLASMA
    plt_VarsName.push_back("I_R(nE)");
#endif
    plt_VarsName.push_back("FunctCall");
    // Extras:
    if (m_plotHeatRelease != 0) {
      plt_VarsName.push_back("HeatRelease");
    }
  }

#ifdef AMREX_USE_EB
  plt_VarsName.push_back("volFrac");
#endif

  for (int ivar = 0; ivar < m_derivePlotVarCount; ++ivar) {
    const PeleLMDeriveRec* rec = derive_lst.get(m_derivePlotVars[ivar]);
    for (int dvar = 0; dvar < rec->numDerive(); ++dvar) {
      plt_VarsName.push_back(rec->variableName(dvar));
    }
  }
#ifdef PELE_USE_SPRAY
  if (SprayParticleContainer::NumDeriveVars() > 0) {
    // We need virtual particles for the lower levels
    setupVirtualParticles(0);
    for (const auto& spray_derive_name :
         SprayParticleContainer::DeriveVarNames()) {
      plt_VarsName.push_back(spray_derive_name);
    }
  }
  if (do_spray_particles && SprayParticleContainer::plot_spray_src) {
    plt_VarsName.push_back("spray_mass_src");
    plt_VarsName.push_back("spray_energy_src");
    AMREX_D_TERM(plt_VarsName.push_back("spray_momentumX_src");
                 , plt_VarsName.push_back("spray_momentumY_src");
                 , plt_VarsName.push_back("spray_momentumZ_src"));
    for (const auto& spray_fuel_name :
         SprayParticleContainer::m_sprayDepNames) {
      plt_VarsName.push_back("spray_" + spray_fuel_name + "_src");
    }
  }
#endif

#ifdef PELE_USE_PLASMA
  if (m_do_extraEFdiags) {
    for (int ivar = 0; ivar < NUM_IONS; ++ivar) {
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        const std::string dir = (idim == 0) ? "X" : ((idim == 1) ? "Y" : "Z");
        plt_VarsName.push_back(
          "DriftFlux_" + names[NUM_SPECIES - NUM_IONS + ivar] + "_" + dir);
      }
    }
  }
#endif

  if (m_do_les && m_plot_les) {
    plt_VarsName.push_back("viscturb");
  }

#if NUM_ODE > 0
  for (int n = 0; n < NUM_ODE; ++n) {
    plt_VarsName.push_back(m_ode_names[n]);
  }
#endif

  // External source terms
  if (m_plot_extSource) {
    for (int ivar = 0; ivar < NVAR; ++ivar) {
      plt_VarsName.push_back("extsource_" + stateVariableName(ivar));
    }
  }

  //----------------------------------------------------------------
  // Fill the plot MultiFabs
  for (int lev = 0; lev <= finest_level; ++lev) {
    int cnt = 0;
    if (m_incompressible != 0) {
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, 0, cnt, AMREX_SPACEDIM, 0);
      cnt += AMREX_SPACEDIM;
    } else {
      // Velocity and density
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, 0, cnt, AMREX_SPACEDIM + 1,
        0);
      cnt += AMREX_SPACEDIM + 1;
      // Species only if requested
      if (m_plotStateSpec != 0) {
        amrex::MultiFab::Copy(
          mf_plt[lev], m_leveldata_new[lev]->state, FIRSTSPEC, cnt, NUM_SPECIES,
          0);
        cnt += NUM_SPECIES;
      }
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, RHOH, cnt, 3, 0);
      cnt += 3;
#ifdef PELE_USE_PLASMA
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, NE, cnt, 2, 0);
      cnt += 2;
#endif
#ifdef PELE_USE_SOOT
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->state, FIRSTSOOT, cnt, NUMSOOTVAR,
        0);
      cnt += NUMSOOTVAR;
#endif
#ifdef PELE_USE_RADIATION
      if (do_rad_solve) {
        amrex::MultiFab::Copy(mf_plt[lev], rad_model->G()[lev], 0, cnt, 1, 0);
        cnt += 1;
        amrex::MultiFab::Copy(
          mf_plt[lev], rad_model->kappa()[lev], 0, cnt, 1, 0);
        cnt += 1;
        amrex::MultiFab::Copy(
          mf_plt[lev], rad_model->emis()[lev], 0, cnt, 1, 0);
        cnt += 1;
      }
#endif
      if (m_has_divu != 0) {
        amrex::MultiFab::Copy(
          mf_plt[lev], m_leveldata_new[lev]->divu, 0, cnt, 1, 0);
        cnt += 1;
      }
    }
    if (m_plot_grad_p != 0) {
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->gp, 0, cnt, AMREX_SPACEDIM, 0);
      cnt += AMREX_SPACEDIM;
    }

    if (m_nAux > 0) {
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldata_new[lev]->auxiliaries, 0, cnt, m_nAux, 0);
      cnt += m_nAux;
    }

    if ((m_do_react != 0) && (m_skipInstantRR == 0) && (m_plot_react != 0)) {
      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldatareact[lev]->I_R, 0, cnt, nCompIR(), 0);
      cnt += nCompIR();

      amrex::MultiFab::Copy(
        mf_plt[lev], m_leveldatareact[lev]->functC, 0, cnt, 1, 0);
      cnt += 1;

      if (m_plotHeatRelease != 0) {
        std::unique_ptr<amrex::MultiFab> mf;
        mf = std::make_unique<amrex::MultiFab>(grids[lev], dmap[lev], 1, 0);
        getHeatRelease(lev, mf.get());
        amrex::MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, 1, 0);
        cnt += 1;
      }
    }

#ifdef AMREX_USE_EB
    amrex::MultiFab::Copy(
      mf_plt[lev], EBFactory(lev).getVolFrac(), 0, cnt, 1, 0);
    cnt += 1;
#endif

    for (int ivar = 0; ivar < m_derivePlotVarCount; ++ivar) {
      std::unique_ptr<amrex::MultiFab> mf;
      mf = derive(m_derivePlotVars[ivar], m_cur_time, lev, 0);
      amrex::MultiFab::Copy(mf_plt[lev], *mf, 0, cnt, mf->nComp(), 0);
      cnt += mf->nComp();
    }
#ifdef PELE_USE_SPRAY
    if (SprayParticleContainer::NumDeriveVars() > 0) {
      const int num_spray_derive = SprayParticleContainer::NumDeriveVars();
      mf_plt[lev].setVal(0., cnt, num_spray_derive);
      SprayPC->computeDerivedVars(mf_plt[lev], lev, cnt);
      if (lev < finest_level) {
        amrex::MultiFab tmp_plt(
          grids[lev], dmap[lev], num_spray_derive, 0, amrex::MFInfo(),
          Factory(lev));
        tmp_plt.setVal(0.);
        VirtPC->computeDerivedVars(tmp_plt, lev, 0);
        amrex::MultiFab::Add(mf_plt[lev], tmp_plt, 0, cnt, num_spray_derive, 0);
      }
      cnt += num_spray_derive;
    }
    if (do_spray_particles && SprayParticleContainer::plot_spray_src) {
      SprayComps scomps = SprayParticleContainer::getSprayComps();
      amrex::MultiFab::Copy(
        mf_plt[lev], *m_spraysource[lev], scomps.rhoSrcIndx, cnt++, 1, 0);
      amrex::MultiFab::Copy(
        mf_plt[lev], *m_spraysource[lev], scomps.engSrcIndx, cnt++, 1, 0);
      amrex::MultiFab::Copy(
        mf_plt[lev], *m_spraysource[lev], scomps.momSrcIndx, cnt,
        AMREX_SPACEDIM, 0);
      cnt += AMREX_SPACEDIM;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        amrex::MultiFab::Copy(
          mf_plt[lev], *m_spraysource[lev], scomps.specSrcIndx + spf, cnt++, 1,
          0);
      }
    }
#endif
#ifdef PELE_USE_PLASMA
    if (m_do_extraEFdiags) {
      amrex::MultiFab::Copy(
        mf_plt[lev], *m_ionsFluxes[lev], 0, cnt, m_ionsFluxes[lev]->nComp(), 0);
      cnt += m_ionsFluxes[lev]->nComp();
    }
#endif

    if (m_do_les && m_plot_les) {
      constexpr amrex::Real fact = 0.5 / AMREX_SPACEDIM;
      auto const& plot_arr = mf_plt[lev].arrays();
      AMREX_D_TERM(
        auto const& mut_arr_x =
          m_leveldata_new[lev]->visc_turb_fc[0].const_arrays();
        , auto const& mut_arr_y =
            m_leveldata_new[lev]->visc_turb_fc[1].const_arrays();
        , auto const& mut_arr_z =
            m_leveldata_new[lev]->visc_turb_fc[2].const_arrays();)
      // interpolate turbulent viscosity from faces to centers
      amrex::ParallelFor(
        mf_plt[lev], [plot_arr, cnt, mut_arr_x, mut_arr_y
#if (AMREX_SPACEDIM == 3)
                      ,
                      mut_arr_z
#endif
      ] AMREX_GPU_DEVICE(int box_no, int i, int j, int k) noexcept {
          plot_arr[box_no](i, j, k, cnt) =
            fact *
            (AMREX_D_TERM(
              mut_arr_x[box_no](i, j, k) + mut_arr_x[box_no](i + 1, j, k),
              +mut_arr_y[box_no](i, j, k) + mut_arr_y[box_no](i, j + 1, k),
              +mut_arr_z[box_no](i, j, k) + mut_arr_z[box_no](i, j, k + 1)));
        });
      amrex::Gpu::streamSynchronize();
      cnt += 1;
    }

#if NUM_ODE > 0
    amrex::MultiFab::Copy(
      mf_plt[lev], m_leveldata_new[lev]->state, FIRSTODE, cnt, NUM_ODE, 0);
    cnt += NUM_ODE;
#endif

    if (m_plot_extSource) {
      amrex::MultiFab::Copy(mf_plt[lev], *m_extSource[lev], 0, cnt, NVAR, 0);
    }

#ifdef AMREX_USE_EB
    if (m_plot_zeroEBcovered != 0) {
      EB_set_covered(mf_plt[lev], 0.0);
    }
#endif
  }

  // No SubCycling, all levels the same step.
  amrex::Vector<int> istep(finest_level + 1, m_nstep);

#ifdef AMREX_USE_HDF5
  if (m_write_hdf5_pltfile) {
    amrex::WriteMultiLevelPlotfileHDF5(
      plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt), plt_VarsName,
      Geom(), m_cur_time, istep, refRatio());
  } else
#endif
  {
    amrex::WriteMultiLevelPlotfile(
      plotfilename, finest_level + 1, GetVecOfConstPtrs(mf_plt), plt_VarsName,
      Geom(), m_cur_time, istep, refRatio());
  }

#ifdef PELE_USE_SPRAY
  if (do_spray_particles) {
    constexpr bool is_spraycheck = false;
    for (int lev = 0; lev <= finest_level; ++lev) {
      SprayPC->SprayParticleIO(lev, is_spraycheck, plotfilename);
      // Remove virtual particles that were made for derived variables
      removeVirtualParticles(lev);
    }
  }
#endif
}

void
PeleLM::WriteHeader(const std::string& name, const bool is_checkpoint) const
{
  if (amrex::ParallelDescriptor::IOProcessor()) {
    std::string HeaderFileName(name + "/Header");
    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
    std::ofstream HeaderFile;

    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

    HeaderFile.open(
      HeaderFileName.c_str(),
      std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);

    if (!HeaderFile.good()) {
      amrex::FileOpenFailed(HeaderFileName);
    }

    HeaderFile.precision(17);
    if (is_checkpoint) {
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
    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      HeaderFile << Geom(0).ProbLo(i) << ' ';
    }
    HeaderFile << '\n';

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      HeaderFile << Geom(0).ProbHi(i) << ' ';
    }
    HeaderFile << '\n';

    // BoxArray
    for (int lev = 0; lev <= finest_level; ++lev) {
      boxArray(lev).writeOn(HeaderFile);
      HeaderFile << '\n';
    }

    // Ambient pressure and typvals
    HeaderFile << m_pNew << "\n";
    for (double typical_value : typical_values) {
      HeaderFile << typical_value << "\n";
    }
  }
}

void
PeleLM::WriteCheckPointFile()
{
  BL_PROFILE("PeleLMeX::WriteCheckPointFile()");

  const std::string& checkpointname =
    amrex::Concatenate(m_check_file, m_nstep, m_ioDigits);

  if (m_verbose != 0) {
    amrex::Print() << "\n Writing checkpoint file: " << checkpointname << "\n";
  }

  //----------------------------------------------------------------
  // Delete checkfiles if present and requested (and have same name)
  if (m_check_overwrite) {
    if (amrex::ParallelContext::IOProcessorSub()) {
      if (amrex::FileExists(checkpointname)) {
        amrex::FileSystem::RemoveAll(checkpointname);
      }
    }
  }

  amrex::VisMF::SetNOutFiles(m_nfiles);

  amrex::PreBuildDirectorHierarchy(
    checkpointname, level_prefix, finest_level + 1, true);

  constexpr bool is_checkpoint = true;
  WriteHeader(checkpointname, is_checkpoint);
  WriteJobInfo(checkpointname);

  for (int lev = 0; lev <= finest_level; ++lev) {
    amrex::VisMF::Write(
      m_leveldata_new[lev]->state,
      amrex::MultiFabFileFullPrefix(
        lev, checkpointname, level_prefix, "state"));

    amrex::VisMF::Write(
      m_leveldata_new[lev]->gp, amrex::MultiFabFileFullPrefix(
                                  lev, checkpointname, level_prefix, "gradp"));

    amrex::VisMF::Write(
      m_leveldata_new[lev]->press,
      amrex::MultiFabFileFullPrefix(lev, checkpointname, level_prefix, "p"));

    if (m_nAux > 0) {
      amrex::VisMF::Write(
        m_leveldata_new[lev]->auxiliaries,
        amrex::MultiFabFileFullPrefix(
          lev, checkpointname, level_prefix, "aux"));
    }

    if (m_incompressible == 0) {
      if (m_has_divu != 0) {
        amrex::VisMF::Write(
          m_leveldata_new[lev]->divu,
          amrex::MultiFabFileFullPrefix(
            lev, checkpointname, level_prefix, "divU"));
      }

      if (m_do_react != 0) {
        amrex::VisMF::Write(
          m_leveldatareact[lev]->I_R,
          amrex::MultiFabFileFullPrefix(
            lev, checkpointname, level_prefix, "I_R"));
      }
    }
  }
#ifdef PELE_USE_SPRAY
  if (do_spray_particles) {
    constexpr bool is_spraycheck = true;
    for (int lev = 0; lev <= finest_level; ++lev) {
      SprayPC->SprayParticleIO(lev, is_spraycheck, checkpointname);
    }
  }
#endif
}

void
PeleLM::ReadCheckPointFile()
{
  BL_PROFILE("PeleLMeX::ReadCheckPointFile()");

  amrex::Print() << "Restarting from checkpoint " << m_restart_chkfile << "\n";

  amrex::Real prob_lo[AMREX_SPACEDIM];
  amrex::Real prob_hi[AMREX_SPACEDIM];

  /***************************************************************************
  ** Load header: set up problem domain (including BoxArray)                 *
  **              allocate PeleLMeX memory (PeleLM::AllocateArrays)            *
  **              (by calling MakeNewLevelFromScratch)                       *
  ****************************************************************************/

  std::string File(m_restart_chkfile + "/Header");

  amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());

  amrex::Vector<char> fileCharPtr;
  amrex::ParallelDescriptor::ReadAndBcastFile(File, fileCharPtr);
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
  finest_level = amrex::min(chk_finest_level, max_level);

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

  // Low coordinates of domain bounding box
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      prob_lo[i++] = std::stod(word);
    }
  }

  // High coordinates of domain bounding box
  std::getline(is, line);
  {
    std::istringstream lis(line);
    int i = 0;
    while (lis >> word) {
      prob_hi[i++] = std::stod(word);
    }
  }

  // Set up problem domain
  amrex::RealBox rb(prob_lo, prob_hi);
  amrex::Geometry::ResetDefaultProbDomain(rb);
  for (int lev = 0; lev <= max_level; ++lev) {
    SetGeometry(
      lev,
      amrex::Geometry(
        Geom(lev).Domain(), rb, Geom(lev).CoordInt(), Geom(lev).isPeriodic()));
  }

  for (int lev = 0; lev <= finest_level; ++lev) {
    // read in level 'lev' BoxArray from Header
    amrex::BoxArray ba;
    ba.readFrom(is);
    GotoNextLine(is);

    // Create distribution mapping
    amrex::DistributionMapping dm{ba, amrex::ParallelDescriptor::NProcs()};
    MakeNewLevelFromScratch(lev, m_cur_time, ba, dm);
  }

  for (int lev = finest_level + 1; lev <= chk_finest_level; ++lev) {
    // read dummy level 'lev' BoxArray if restarting with reduced levels
    amrex::BoxArray ba;
    ba.readFrom(is);
  }

  // deal with typval and P_amb
  is >> m_pNew;
  GotoNextLine(is);
  m_pOld = m_pNew;
  for (double& typical_value : typical_values) {
    is >> typical_value;
    GotoNextLine(is);
  }

  /***************************************************************************
   * Load fluid data                                                         *
   ***************************************************************************/

  // Load the field data
  for (int lev = 0; lev <= finest_level; ++lev) {
#ifdef PELE_USE_PLASMA
    if (!m_restart_nonEF) {
      amrex::VisMF::Read(
        m_leveldata_new[lev]->state,
        amrex::MultiFabFileFullPrefix(
          lev, m_restart_chkfile, level_prefix, "state"));
    } else {
      // The chk state is 2 component shorter since phiV and nE aren't in it
      amrex::MultiFab stateTemp(grids[lev], dmap[lev], NVAR - 2, m_nGrowState);
      amrex::VisMF::Read(
        stateTemp, amrex::MultiFabFileFullPrefix(
                     lev, m_restart_chkfile, level_prefix, "state"));
      amrex::MultiFab::Copy(
        m_leveldata_new[lev]->state, stateTemp, 0, 0, NVAR - 2, m_nGrowState);
    }
#else
    amrex::VisMF::Read(
      m_leveldata_new[lev]->state,
      amrex::MultiFabFileFullPrefix(
        lev, m_restart_chkfile, level_prefix, "state"));
#endif

    amrex::VisMF::Read(
      m_leveldata_new[lev]->gp,
      amrex::MultiFabFileFullPrefix(
        lev, m_restart_chkfile, level_prefix, "gradp"));

    amrex::VisMF::Read(
      m_leveldata_new[lev]->press,
      amrex::MultiFabFileFullPrefix(lev, m_restart_chkfile, level_prefix, "p"));
    if (m_nAux > 0) {
      amrex::VisMF::Read(
        m_leveldata_new[lev]->auxiliaries,
        amrex::MultiFabFileFullPrefix(
          lev, m_restart_chkfile, level_prefix, "aux"));
    }

    if (m_incompressible == 0) {
      if (m_has_divu != 0) {
        amrex::VisMF::Read(
          m_leveldata_new[lev]->divu,
          amrex::MultiFabFileFullPrefix(
            lev, m_restart_chkfile, level_prefix, "divU"));
      }

#ifdef PELE_USE_PLASMA
      if (!m_restart_nonEF) {
        if (m_do_react) {
          amrex::VisMF::Read(
            m_leveldatareact[lev]->I_R,
            amrex::MultiFabFileFullPrefix(
              lev, m_restart_chkfile, level_prefix, "I_R"));
        }
      } else {
        // I_R for non-EF simulation is one component shorted, need to account
        // for that.
        if (m_do_react) {
          amrex::MultiFab I_Rtemp(grids[lev], dmap[lev], NUM_SPECIES, 0);
          amrex::VisMF::Read(
            I_Rtemp, amrex::MultiFabFileFullPrefix(
                       lev, m_restart_chkfile, level_prefix, "I_R"));
          amrex::MultiFab::Copy(
            m_leveldatareact[lev]->I_R, I_Rtemp, 0, 0, NUM_SPECIES, 0);
        }

        // Initialize nE & phiV
        m_leveldata_new[lev]->state.setVal(0.0, NE, 2, m_nGrowState);
      }
#else
      if (m_do_react != 0) {
        amrex::VisMF::Read(
          m_leveldatareact[lev]->I_R,
          amrex::MultiFabFileFullPrefix(
            lev, m_restart_chkfile, level_prefix, "I_R"));
      }
#endif
    }
  }
  if (m_verbose != 0) {
    amrex::Print() << "Restart complete \n";
  }
}

void
PeleLM::initLevelDataFromPlt(int a_lev, const std::string& a_dataPltFile)
{
  if (m_incompressible != 0) {
    amrex::Abort(
      " initializing data from a pltfile only available for low-Mach "
      "simulations");
  }
  if (m_nAux > 0) {
    amrex::Warning(
      " restarting from plotfile with auxiliaries not currently "
      "implemented, and will not be captured");
  }
  if (m_verbose > 0) {
    amrex::Print() << " initData on level " << a_lev << " from pltfile "
                   << a_dataPltFile << "\n";
    if (pltfileSource == "LM") {
      amrex::Print() << " Assuming pltfile was generated in LM/LMeX \n";
    } else if (pltfileSource == "C") {
      amrex::Print() << " Assuming pltfile was generated in PeleC \n";
    }
  }

  // Use PelePhysics PltFileManager
  pele::physics::pltfilemanager::PltFileManager pltData(a_dataPltFile);
  amrex::Vector<std::string> plt_vars = pltData.getVariableList();
  if (m_do_reset_time == 0) {
    m_cur_time = pltData.getTime();
    m_nstep = pltData.getNsteps();
  }

  // Find required data in pltfile
  amrex::Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names, &(eos_parms.host_parm()));
  int idT = -1, idV = -1, idY = -1, nSpecPlt = 0;
#ifdef PELE_USE_PLASMA
  int inE = -1, iPhiV = -1;
#endif
#ifdef PELE_USE_SOOT
  int inSoot = -1;
#endif
  for (int i = 0; i < plt_vars.size(); ++i) {
    std::string firstChars = plt_vars[i].substr(0, 2);

    if (pltfileSource == "LM") {
      if (plt_vars[i] == "temp") {
        idT = i;
      }
    } else if (pltfileSource == "C") {
      if (plt_vars[i] == "Temp") {
        idT = i;
      }
    }
    if (plt_vars[i] == "x_velocity") {
      idV = i;
    }
    if (firstChars == "Y(" && idY < 0) { // species might not be ordered in the
                                         // order of the current mech.
      idY = i;
    }
    if (firstChars == "Y(") {
      nSpecPlt += 1;
    }
#ifdef PELE_USE_PLASMA
    if (plt_vars[i] == "nE")
      inE = i;
    if (plt_vars[i] == "phiV")
      iPhiV = i;
#endif
#ifdef PELE_USE_SOOT
    if (plt_vars[i] == "soot_N") {
      inSoot = i;
    }
#endif
  }
  if (idY < 0) {
    amrex::Abort("Couldn't find species mass fractions in pltfile");
  } else if (idT < 0) {
    amrex::Abort("Couldn't find temperature in pltfile");
  }
  if (m_verbose > 0) {
    amrex::Print() << " " << nSpecPlt
                   << " species found in pltfile, starting with "
                   << plt_vars[idY] << "\n";
  }

  // Get level data
  auto* ldata_p = getLevelDataPtr(a_lev, AmrNewTime);

  // Velocity
  pltData.fillPatchFromPlt(
    a_lev, geom[a_lev], idV, VELX, AMREX_SPACEDIM, ldata_p->state);

  // Temperature
  pltData.fillPatchFromPlt(a_lev, geom[a_lev], idT, TEMP, 1, ldata_p->state);

  // Species
  // Hold the species in temporary MF before copying to level data
  // in case the number of species differs.
  amrex::MultiFab speciesPlt(grids[a_lev], dmap[a_lev], nSpecPlt, 0);
  pltData.fillPatchFromPlt(a_lev, geom[a_lev], idY, 0, nSpecPlt, speciesPlt);
  for (int i = 0; i < NUM_SPECIES; ++i) {
    std::string specName = "Y(" + spec_names[i] + ")";
    std::string specString = specName;
    if (!m_initDataPlt_specname_map.empty()) {
      specString = m_initDataPlt_specname_map[i];
    }
    int foundSpec = 0;
    for (int iplt = 0; iplt < nSpecPlt; ++iplt) {
      if (specString == plt_vars[idY + iplt]) {
        amrex::MultiFab::Copy(
          ldata_p->state, speciesPlt, iplt, FIRSTSPEC + i, 1, 0);
        foundSpec = 1;
        if (m_verbose > 0) {
          amrex::Print() << "Loading species " << specName
                         << " from plotfile species " << specString << "\n";
        }
      }
    }
    if (foundSpec == 0) {
      for (int iplt = 0; iplt < plt_vars.size(); ++iplt) {
        if (specString == plt_vars[iplt]) {
          foundSpec = 1;
          if (m_verbose > 0) {
            amrex::Print() << "Loading species " << specName
                           << " from plotfile entry " << specString << "\n";
          }
          pltData.fillPatchFromPlt(
            a_lev, geom[a_lev], iplt, FIRSTSPEC + i, 1, ldata_p->state);
        }
      }
    }
    if (foundSpec == 0) {
      ldata_p->state.setVal(0.0, FIRSTSPEC + i, 1);
      if (m_verbose > 0) {
        amrex::Print() << "For species " << specName << " entry " << specString
                       << " not found in plot file, setting to 0\n";
      }
    }
  }

  // Converting units when pltfile is coming from PeleC solution
  if (pltfileSource == "C") {
    if (m_verbose > 0) {
      amrex::Print() << " Converting CGS to MKS units... \n";
    }
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(ldata_p->state, amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      auto const& vel_arr = ldata_p->state.array(mfi, VELX);
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
            amrex::Real vel_mks = vel_arr(i, j, k, n) * 0.01;
            vel_arr(i, j, k, n) = vel_mks;
          }
        });
    }
  }

#ifdef PELE_USE_PLASMA
  // nE
  pltData.fillPatchFromPlt(a_lev, geom[a_lev], inE, NE, 1, ldata_p->state);
  // phiV
  pltData.fillPatchFromPlt(a_lev, geom[a_lev], iPhiV, PHIV, 1, ldata_p->state);
#endif
#ifdef PELE_USE_SOOT
  if (do_soot_solve) {
    if (inSoot >= 0) {
      pltData.fillPatchFromPlt(
        a_lev, geom[a_lev], inSoot, FIRSTSOOT, NUMSOOTVAR, ldata_p->state);
      if (pltfileSource == "C") {
        SootConst sc;
        amrex::Real* momV = sc.MomOrderV.data();
        amrex::Real* momS = sc.MomOrderS.data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
        for (amrex::MFIter mfi(ldata_p->state, amrex::TilingIfNotGPU());
             mfi.isValid(); ++mfi) {
          const amrex::Box& bx = mfi.tilebox();
          auto const& soot_arr = ldata_p->state.array(mfi, FIRSTSOOT);
          amrex::ParallelFor(
            bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              for (int n = 0; n < NUM_SOOT_MOMENTS; n++) {
                amrex::Real soot_exp = 3. - (3. * momV[n] + 2. * momS[n]);
                soot_arr(i, j, k, n) *= std::pow(100., soot_exp);
              }
              soot_arr(i, j, k, NUMSOOTVAR - 1) *= 1.E6;
            });
        }
      }
    } else {
      SootData* const sd = soot_model->getSootData();
      amrex::Real moments[NUM_SOOT_MOMENTS + 1];
      sd->initialSmallMomVals(moments);
      for (int mom = 0; mom < NUM_SOOT_MOMENTS + 1; ++mom) {
        ldata_p->state.setVal(moments[mom], FIRSTSOOT + mom, 1);
      }
    }
  }
#endif
  // Pressure and pressure gradients to zero
  ldata_p->press.setVal(0.0);
  ldata_p->gp.setVal(0.0);

  ProbParm const* lprobparm = prob_parm_d;
  auto const* leosparm = eos_parms.device_parm();

  // If m_do_patch_flow_variables is set as true, call user-defined function to
  // patch flow variables
  if (m_do_patch_flow_variables) {
    ProblemSpecificFunctions::patchFlowVariables(
      geom[a_lev], *lprobparm, ldata_p->state);
  }

  // Enforce rho and rhoH consistent with temperature and mixture
  // The above handles species mapping (to some extent), but nothing enforce
  // sum of Ys = 1 -> use N2 in the following if N2 is present
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(ldata_p->state, amrex::TilingIfNotGPU());
       mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    auto const& rho_arr = ldata_p->state.array(mfi, DENSITY);
    auto const& rhoY_arr = ldata_p->state.array(mfi, FIRSTSPEC);
    auto const& rhoH_arr = ldata_p->state.array(mfi, RHOH);
    auto const& temp_arr = ldata_p->state.array(mfi, TEMP);
    const auto* eosparm = leosparm;
    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      auto eos = pele::physics::PhysicsType::eos(eosparm);
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      amrex::Real sumYs = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = rhoY_arr(i, j, k, n);
#ifdef N2_ID
        if (n != N2_ID) {
          sumYs += massfrac[n];
        }
#endif
      }
#ifdef N2_ID
      massfrac[N2_ID] = 1.0 - sumYs;
#endif

      // Get density
      amrex::Real P_cgs = lprobparm->P_mean * 10.0;
      amrex::Real rho_cgs = 0.0;
      eos.PYT2R(P_cgs, massfrac, temp_arr(i, j, k), rho_cgs);
      rho_arr(i, j, k) = rho_cgs * 1.0e3;

      // Get enthalpy
      amrex::Real h_cgs = 0.0;
      eos.TY2H(temp_arr(i, j, k), massfrac, h_cgs);
      rhoH_arr(i, j, k) = h_cgs * 1.0e-4 * rho_arr(i, j, k);

      // Fill rhoYs
      for (int n = 0; n < NUM_SPECIES; n++) {
        rhoY_arr(i, j, k, n) = massfrac[n] * rho_arr(i, j, k);
      }
    });
  }

  // Initialize thermodynamic pressure
  setThermoPress(a_lev, AmrNewTime);
  if (m_has_divu != 0) {
    ldata_p->divu.setVal(0.0);
  }
}

void
PeleLM::addLevelVelocityDataFromPlt(int a_lev, const std::string& a_velPltFile)
{
  if (m_verbose > 0) {
    amrex::Print() << " init velocity data on level " << a_lev
                   << " from pltfile " << a_velPltFile << "\n";
  }

  // Use PelePhysics PltFileManager
  pele::physics::pltfilemanager::PltFileManager pltData(a_velPltFile);
  amrex::Vector<std::string> plt_vars = pltData.getVariableList();

  // do some compatibility checks
  if (pltData.getNlev() < a_lev) {
    amrex::Abort("USE_VELOCITY: not enough levels in plotfile");
  }
  if (pltData.getGeom(a_lev).Domain() != geom[a_lev].Domain()) {
    amrex::Abort("USE_VELOCITY: problem domains do not match");
  }

  // find velocity in the plotfile
  int idXvel = -1;
  for (int i = 0; i < plt_vars.size(); ++i) {
    if (plt_vars[i] == "x_velocity") {
      idXvel = i;
    }
  }
  if (idXvel == -1) {
    amrex::Abort(
      "Could not find velocity fields in supplied velocity_plotfile");
  }

  // Get level data
  auto* ldata_p = getLevelDataPtr(a_lev, AmrNewTime);

  // load data from plot file
  amrex::BoxArray tmpVelBA(ldata_p->state.boxArray());
  amrex::DistributionMapping tmpVelDM(tmpVelBA);
  int nGrow0(0), sComp0(0);
  amrex::MultiFab tmpVel(tmpVelBA, tmpVelDM, AMREX_SPACEDIM, nGrow0);
  pltData.fillPatchFromPlt(
    a_lev, geom[a_lev], idXvel, sComp0, AMREX_SPACEDIM, tmpVel);
  // scale the velocity
  tmpVel.mult(m_velocity_plotfile_scale);
  amrex::MultiFab::Add(ldata_p->state, tmpVel, 0, VELX, AMREX_SPACEDIM, 0);
}

void
PeleLM::WriteJobInfo(const std::string& path) const
{
  if (amrex::ParallelDescriptor::IOProcessor()) {
    // job_info file with details about the run
    std::ofstream jobInfoFile;
    std::string FullPathJobInfoFile = path;

    FullPathJobInfoFile += "/PeleLMeX_job_info";
    jobInfoFile.open(FullPathJobInfoFile.c_str(), std::ios::out);

    // job information
    jobInfoFile << PrettyLine;
    jobInfoFile << " PeleLMeX Job Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "number of MPI processes: "
                << amrex::ParallelDescriptor::NProcs() << "\n";
#ifdef AMREX_USE_OMP
    jobInfoFile << "number of threads:       " << omp_get_max_threads() << "\n";
#endif

    jobInfoFile << "\n\n";

    // build information
    jobInfoFile << PrettyLine;
    jobInfoFile << " Build Information\n";
    jobInfoFile << PrettyLine;

    jobInfoFile << "build date:    " << amrex::buildInfoGetBuildDate() << "\n";
    jobInfoFile << "build machine: " << amrex::buildInfoGetBuildMachine()
                << "\n";
    jobInfoFile << "build dir:     " << amrex::buildInfoGetBuildDir() << "\n";
    jobInfoFile << "AMReX dir:     " << amrex::buildInfoGetAMReXDir() << "\n";

    jobInfoFile << "\n";

    jobInfoFile << "COMP:          " << amrex::buildInfoGetComp() << "\n";
    jobInfoFile << "COMP version:  " << amrex::buildInfoGetCompVersion()
                << "\n";
    jobInfoFile << "C++ compiler:  " << amrex::buildInfoGetCXXName() << "\n";
    jobInfoFile << "C++ flags:     " << amrex::buildInfoGetCXXFlags() << "\n";

    jobInfoFile << "\n";

    const char* githash1 = amrex::buildInfoGetGitHash(1);
    const char* githash2 = amrex::buildInfoGetGitHash(2);
    const char* githash3 = amrex::buildInfoGetGitHash(3);
    const char* githash4 = amrex::buildInfoGetGitHash(4);

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

    for (int lev = 0; lev <= finest_level; ++lev) {
      jobInfoFile << " level: " << lev << "\n";
      jobInfoFile << "   number of boxes = " << grids[lev].size() << "\n";
      jobInfoFile << "   maximum zones   = ";
      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        jobInfoFile << geom[lev].Domain().length(idim) << " ";
      }
      jobInfoFile << "\n\n";
    }

    jobInfoFile << "\n\n";

    // runtime parameters
    jobInfoFile << PrettyLine;
    jobInfoFile << " Inputs File Parameters\n";
    jobInfoFile << PrettyLine;

    amrex::ParmParse::dumpTable(jobInfoFile, true);

    jobInfoFile.close();
  }
}
