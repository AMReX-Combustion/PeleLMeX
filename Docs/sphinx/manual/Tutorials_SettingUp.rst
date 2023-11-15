Setting-up your environment
---------------------------

Getting a functioning environment in which to compile and run `PeleLMeX` is the first step of this tutorial.
Follow the steps listed below to get to this point:

#. The first step is to get `PeleLMeX` and its dependencies. To do so, use a recursive *git clone*: ::

    git clone --recursive https://github.com/AMReX-Combustion/PeleLMeX.git

#. Move into the Exec folder containing the ``EB_BackwardStepFlame``. To do so: ::

    cd PeleLMeX/Exec/RegTests/EB_BackwardStepFlame

Note that the makefile system is set up such that default paths are automatically set to the
submodules obtained with the recursive *git clone*, however the user can set their own dependencies
in the `GNUmakefile` by updating the top-most lines as follows: ::

       PELELMEX_HOME     = <path_to_PeleLMeX>
       AMREX_HOME        = <path_to_MyAMReX>
       AMREX_HYDRO_HOME  = <path_to_MyAMReXHydro>
       PELE_PHYSICS_HOME = <path_to_MyPelePhysics>
       SUNDIALS_HOME     = <path_to_MySUNDIALS>

or directly through shell environment variables (using *bash* for instance): ::

       export PELELMEX_HOME=<path_to_PeleLMeX>
       export AMREX_HOME=<path_to_MyAMReX>
       export AMREX_HYDRO_HOME=<path_to_MyAMReXHydro>
       export PELE_PHYSICS_HOME=<path_to_MyPelePhysics>
       export SUNDIALS_HOME=<path_to_MySUNDIALS>

Note that using the first option will overwrite any
environment variables you might have previously defined when using this `GNUmakefile`.

You're good to go !
