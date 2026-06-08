#include <PeleLMeX.H>

#ifdef AMREX_USE_ASCENT

#include <AMReX_Conduit_Blueprint.H>
#include <ascent.hpp>

void
PeleLM::doInSituViz()
{
  BL_PROFILE("PeleLMeX::doInSituViz()");

  if (m_ascent.plot_int <= 0) {
    return;
  }
  if (m_nstep % m_ascent.plot_int != 0) {
    return;
  }

  auto dPlotFileTime0 = amrex::second();

  averageDownState(AmrNewTime);

  amrex::Vector<std::string> speciesNames;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    speciesNames, &(eos_parms.host_parm()));

  amrex::Vector<std::string> plt_var_names;
  AMREX_D_TERM(plt_var_names.push_back("x_velocity");
               , plt_var_names.push_back("y_velocity");
               , plt_var_names.push_back("z_velocity"));
  plt_var_names.push_back("density");
  for (int n = 0; n < NUM_SPECIES; ++n) {
    plt_var_names.push_back("rho.Y(" + speciesNames[n] + ")");
  }
  plt_var_names.push_back("rhoh");
  plt_var_names.push_back("temp");
  plt_var_names.push_back("RhoRT");

  // Base state: velocity + density + species + rhoh + temp + RhoRT
  const int ncomp = RHORT + 1;
  AMREX_ASSERT(ncomp == static_cast<int>(plt_var_names.size()));

  const int nlevels = finest_level + 1;
  amrex::Vector<amrex::MultiFab> plotMFs;
  plotMFs.reserve(nlevels);
  for (int lev = 0; lev <= finest_level; ++lev) {
    plotMFs.emplace_back(
      grids[lev], dmap[lev], ncomp, 0, amrex::MFInfo(), Factory(lev));
    amrex::MultiFab::Copy(
      plotMFs[lev], m_leveldata_new[lev]->state, 0, 0, ncomp, 0);
  }

  amrex::Vector<const amrex::MultiFab*> plotMFs_constvec;
  plotMFs_constvec.reserve(nlevels);
  for (int lev = 0; lev < nlevels; ++lev) {
    plotMFs_constvec.push_back(
      static_cast<const amrex::MultiFab*>(&plotMFs[lev]));
  }

  // PeleLMeX is non-subcycling: all levels share m_nstep
  amrex::Vector<int> istep(nlevels, m_nstep);

  conduit::Node bp_mesh;
  amrex::MultiLevelToBlueprint(
    nlevels, plotMFs_constvec, plt_var_names, Geom(), m_cur_time, istep,
    refRatio(), bp_mesh);

  ascent::Ascent ascent;
  conduit::Node open_opts;
#ifdef AMREX_USE_MPI
  open_opts["mpi_comm"] =
    MPI_Comm_c2f(amrex::ParallelDescriptor::Communicator());
#endif
  ascent.open(open_opts);

  conduit::Node verify_info;
  if (!conduit::blueprint::mesh::verify(bp_mesh, verify_info)) {
    ASCENT_INFO("Error: Mesh Blueprint Verify Failed!");
    verify_info.print();
  }

  conduit::Node actions;
  ascent.publish(bp_mesh);
  ascent.execute(actions);
  ascent.close();

  if (m_verbose > 0) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    auto dPlotFileTime = amrex::second() - dPlotFileTime0;
    amrex::ParallelDescriptor::ReduceRealMax(dPlotFileTime, IOProc);
    amrex::Print() << "Ascent write time = " << dPlotFileTime << "  seconds"
                   << std::endl;
  }
}

#endif // AMREX_USE_ASCENT
