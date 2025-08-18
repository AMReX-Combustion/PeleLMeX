#include <PeleLMeX.H>

// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
#include <AMReX_Sundials.H>
#endif

int
main(int argc, char* argv[])
{
  if (argc >= 2) {
    for (auto i = 1; i < argc; ++i) {
      if (std::string(argv[i]) == "--describe") {
        writeBuildInfo();
        return 0;
      }
    }
  }

  amrex::Initialize(argc, argv);

// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
  amrex::sundials::Initialize(amrex::OpenMP::get_max_threads());
#endif

  if (argc <= 1) {
    amrex::Abort("Error: no inputs file provided on command line.");
  }

  BL_PROFILE_VAR("PeleLMeX::main()", main);

  const amrex::Real strt_total = amrex::ParallelDescriptor::second();

  {
    PeleLM pelelmex;
    pelelmex.Setup();
    pelelmex.Init();

    if (pelelmex.runMode() == "normal") {
      pelelmex.Evolve();
    } else if (pelelmex.runMode() == "evaluate") {
      pelelmex.Evaluate();
    } else {
      amrex::Abort(
        " Wrong peleLM.run_mode ! It can only be 'normal' (D) or 'evaluate'");
    }

    amrex::Real end_total = amrex::ParallelDescriptor::second() - strt_total;

    amrex::ParallelDescriptor::ReduceRealMax(
      end_total, amrex::ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "\nTotal Time: " << end_total << '\n';
  }

  BL_PROFILE_VAR_STOP(main);

// Defined and finalized when in gnumake, but not defined in cmake and
// finalization done manually
#ifndef AMREX_USE_SUNDIALS
  amrex::sundials::Finalize();
#endif

  amrex::Finalize();
}
