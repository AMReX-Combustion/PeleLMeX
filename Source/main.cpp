#include <PeleLM.H>

using namespace amrex;

int main(int argc, char* argv[]) {

   // check to see if it contains --describe
   if (argc >= 2) {
      for (auto i = 1; i < argc; i++) {
         if (std::string(argv[i]) == "--describe") {
            writeBuildInfo();
            return 0;
         }
      }
   }

   // in AMReX.cpp
   Initialize(argc, argv);

   // Refuse to continue if we did not provide an inputs file.
   if (argc <= 1) {
       Abort("Error: no inputs file provided on command line.");
   }

   // timer for profiling
   BL_PROFILE_VAR("PeleLM::main()", main);

   // wallclock time
   const Real strt_total = ParallelDescriptor::second();

   {
       // declare an PeleLM object to manage multilevel data
       PeleLM pelelm;

       // Description Setup
       pelelm.Setup();

       // Description Init
       pelelm.Init();

       // Switch between Evolve and UnitTest mode
       if ( pelelm.runMode() == "normal" ) {

          // Advance solution to final time
          pelelm.Evolve();

       } else if ( pelelm.runMode() == "evaluate" ) {

          //
          pelelm.Evaluate();

       } else {
          Abort(" Wrong peleLM.run_mode ! It can only be 'normal' (D) or 'evaluate'");
       }

       // wallclock time
       Real end_total = ParallelDescriptor::second() - strt_total;

       // print wallclock time
       ParallelDescriptor::ReduceRealMax( end_total, ParallelDescriptor::IOProcessorNumber());
       Print() << "\nTotal Time: " << end_total << '\n';
   }

   // destroy timer for profiling
   BL_PROFILE_VAR_STOP(main);

   // in AMReX.cpp
   Finalize();
}
