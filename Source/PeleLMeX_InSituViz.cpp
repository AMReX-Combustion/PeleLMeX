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

  // Assemble the same plot MultiFab and variable list used by WritePlotFile().
  // This ensures Ascent publishes all user-requested fields (species flag,
  // plasma, soot, radiation, reactions, derived vars, spray, LES viscosity,
  // ODE quantities, EB volume fraction, external sources) rather than only
  // the hardcoded base state.  averageDown is performed inside constructPlotMF.
  // constructPlotMF also calls setupVirtualParticles(0) when PELE_USE_SPRAY is
  // active; we are responsible for the matching removeVirtualParticles cleanup
  // since there is no SprayParticleIO step here (unlike WritePlotFile).
  amrex::Vector<amrex::MultiFab> plotMFs;
  amrex::Vector<std::string> plt_var_names;
  constructPlotMF(plotMFs, plt_var_names);

#ifdef PELE_USE_SPRAY
  if (do_spray_particles) {
    for (int lev = 0; lev <= finest_level; ++lev) {
      removeVirtualParticles(lev);
    }
  }
#endif

  // PeleLMeX is non-subcycling: all levels share m_nstep
  amrex::Vector<int> istep(finest_level + 1, m_nstep);

  conduit::Node bp_mesh;
  amrex::MultiLevelToBlueprint(
    finest_level + 1, GetVecOfConstPtrs(plotMFs), plt_var_names, Geom(),
    m_cur_time, istep, refRatio(), bp_mesh);

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
