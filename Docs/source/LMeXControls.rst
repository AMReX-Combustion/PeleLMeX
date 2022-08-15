PeleLMeX controls
=================

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.
This file needs to specified along with the executable as an `argv` option, for example:

::

    mpirun -np 64 ./PeleLMeX2d.xxx.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inputs file, e.g.:

::

    mpirun -np 64 ./PeleLMeX2d.xxx.ex inputs amr.max_level=2

The available options are divided into groups: those that control primarily AMReX are prefaced with `amr.`, those that are specific to the PeleLM are prefaced by `peleLM.`, while those corresponding to the various pieces of the algorithm are prefaced with specific keys, such that `diffusion`, `nodal_proj`, ... as described below.

Computational domain definition
-------------------------------

::

    #--------------------GEOMETRY DEFINITION-----------------------
    geometry.is_periodic = 1 1 0              # For each dir, 0: non-perio, 1: periodic
    geometry.coord_sys   = 0                  # 0 => cart, 1 => RZ
    geometry.prob_lo     = 0.0   0.0   0.0    # x_lo y_lo (z_lo)
    geometry.prob_hi     = 0.016 0.016 0.032  # x_hi y_hi (z_hi)

    #-----------------------BNDY CONDITIONS------------------------
    # Interior, Inflow, Outflow, Symmetry, SlipWallAdiab, NoSlipWallAdiab
    # SlipWallIsotherm, NoSlipWallIsotherm
    # Periodic direction must be set as Interior
    peleLM.lo_bc = Interior Interior Inflow
    peleLM.hi_bc = Interior Interior Inflow

AMR parameters
--------------

::

    #-------------------------AMR CONTROL--------------------------
    amr.n_cell          = 64 64 128        # Number of cells on Level 0 in each direction   
    amr.v               = 1                # [OPT, DEF=0] AMR verbose
    amr.max_level       = 1                # maximum level number allowed
    amr.ref_ratio       = 2 2 2 2          # refinement ratio, one per refinement level
    amr.regrid_int      = 5                # how often to regrid
    amr.n_error_buf     = 1 1 2 2          # number of buffer cells in error est, one per refinement level
    amr.grid_eff        = 0.7              # what constitutes an efficient grid
    amr.blocking_factor = 16               # block factor in grid generation (min box size)
    amr.max_grid_size   = 64               # max box size

Time stepping parameters
------------------------

::

    #-------------------------TIME STEPPING------------------------
    amr.max_step      = 20                 # Maximum number of steps
    amr.stop_time     = 0.001              # Maximum simulation time [s]
    amr.cfl           = 0.5                # [OPT, DEF=0.7] CFL for advection-controlled dt estimate
    amr.fixed_dt      = 1e-6               # [OPT] optional fixed dt (override CFL condition)
    amr.min_dt        = 1e-11              # [OPT, DEF=1e-12] small time step size limit triggering simulation termination
    amr.init_dt       = 1e-6               # [OPT] optional initial dt (override CFL condition upon initialization)
    amr.dt_shrink     = 0.0001             # [OPT, DEF=1.0] dt factor upon initialization
    amr.dt_change_max = 1.1                # [OPT, DEF=1.1] maximum dt change between consecutive steps

Note that either a `max_step` or a `stop_time` is required, and if both are specified, the first stopping criteria
encountered will lead to termination of the simulation.

IO parameters
-------------

::

    #--------------------------IO CONTROL--------------------------
    amr.plot_int         = 20              # [OPT, DEF=-1] Frequency (as step #) for writting plot file
    amr.plot_per         = 002             # [OPT, DEF=-1] Period (time in s) for writting plot file
    amr.plot_per_exact   = 1               # [OPT, DEF=0] Flag to enforce exactly plt_per by shortening dt 
    amr.plot_file        = "plt_"          # [OPT, DEF="plt_"] Plot file prefix
    amr.check_int        = 100             # [OPT, DEF=-1] Frequency (as step #) for writting checkpoint file
    amr.check_file       = "chk"           # [OPT, DEF="chk"] Checkpoint file prefix
    amr.file_stepDigits  = 6               # [OPT, DEF=5] Number of digits when adding nsteps to plt and chk names
    amr.derive_plot_vars = avg_pressure ...# [OPT, DEF=""] List of derived variable included in the plot files
    amr.plot_speciesState = 0              # [OPT, DEF=0] Force adding state rhoYs to the plot files

    amr.restart          = chk00100        # [OPT, DEF=""] Checkpoint from which to restart the simulation
    amr.initDataPlt      = plt01000        # [OPT, DEF=""] Provide a plotfile from which to extract initial data

PeleLMeX derived variables
--------------------------

The following list of derived variables are available in PeleLMeX:

.. list-table:: PeleLMeX derived variables
    :widths: 25 25 100
    :header-rows: 1

    * - Key
      - Size (nComp)
      - Brief
    * - `mass_fractions`
      - NUM_SPECIES
      - Species mass fractions
    * - `mole_fractions`
      - NUM_SPECIES
      - Species mole fractions
    * - `diffcoeffs`
      - NUM_SPECIES
      - Species mixture-averaged diffusion coefficients
    * - `lambda`
      - 1
      - Thermal diffusivity
    * - `viscosity`
      - 1
      - Mixture viscosity
    * - `mixture_fraction`
      - 1
      - Mixture fraction based on Bilger's element formulation
    * - `progress_variable`
      - 1
      - Progress variable based on a linear combination of Ys, T
    * - `avg_pressure`
      - 1
      - Cell-averaged pressure (from the node-centered pressure)
    * - `mag_vort`
      - 1
      - Vorticity (2D) or vorticity magnitude (3D)
    * - `kinetic_energy`
      - 1
      - Kinetic energy: 0.5 * rho * (u^2+v^2+w^2)
    * - `enstrophy`
      - 1
      - enstrophy: 0.5 * rho * (\omega_x^2+\omega_y^2+\omega_z^2)
    * - `HeatRelease`
      - 1
      - Heat release rate from chem. reactions

Note that `mixture_fraction` and `progress_variable` requires additional inputs from the users as described below.

PeleLMeX algorithm
------------------

::

    #-----------------------PELELMEX CONTROL-----------------------
    peleLM.v = 1                           # [OPT, DEF=0] Verbose
    peleLM.run_mode = normal               # [OPT, DEF=normal] Switch between time-advance mode (normal) or UnitTest (evaluate)
    peleLM.use_wbar = 1                    # [OPT, DEF=1] Enable Wbar correction in diffusion fluxes
    peleLM.sdc_iterMax = 2                 # [OPT, DEF=1] Number of SDC iterations
    peleLM.num_init_iter = 2               # [OPT, DEF=3] Number of iterations to get initial pressure
    peleLM.num_divu_iter = 1               # [OPT, DEF=1] Number of divU iterations to get initial dt estimate
    peleLM.do_init_proj = 1                # [OPT, DEF=1] Control over initial projection
    peleLM.advection_scheme = Godunov_BDS  # [OPT, DEF=Godunov_PLM] Advection scheme: Godunov_PLM, Godunov_PPM or Godunov_BDS
    peleLM.incompressible = 0              # [OPT, DEF=0] Enable to run fully incompressible, scalar advance is bypassed
    peleLM.m_rho = 1.17                    # [OPT, DEF=-1] If incompressible, density value [MKS]
    peleLM.m_mu = 1.8e-5                   # [OPT, DEF=-1] If incompressible, kinematic visc. value [MKS]
    peleLM.gravity = 0.0 0.0 -9.81         # [OPT, DEF=Vec{0.0}] Gravity vector [MKS]
    peleLM.gradP0 = 0.0 0.0 10.0           # [OPT, DEF=Vec{0.0}] Average background pressure gradient [Pa/m]
    peleLM.closed_chamber = 0              # [OPT] Override the automatic detection of closed chamber (based on Outflow(s))
    peleLM.floor_species = 0               # [OPT, DEF=0] Crudely enforce mass fraction positivity
    peleLM.deltaT_verbose = 0              # [OPT, DEF=0] Verbose of the deltaT iterative solve algorithm
    peleLM.deltaT_iterMax = 5              # [OPT, DEF=10] Maximum number of deltaT iterations
    peleLM.deltaT_tol = 1e-10              # [OPT, DEF=1.e-10] Tolerance of the deltaT solve
    peleLM.evaluate_vars =...              # [OPT, DEF=""] In evaluate mode, list unitTest: diffTerm, divU, instRR, transportCC

Chemistry integrator
--------------------

::

    #-----------------------CHEMISTRY CONTROL----------------------
    peleLM.chem_integrator   = "ReactorCvode"   # Chemistry integrator, from PelePhysics available list
    peleLM.use_typ_vals_chem = 1                # [OPT, DEF=1] Use Typical values to scale components in the reactors
    peleLM.typical_values_reset_int = 5         # [OPT, DEF=10] Frequency at which the typical values are updated
    ode.rtol = 1.0e-6                           # [OPT, DEF=1e-10] Relative tolerance of the chem. reactor
    ode.atol = 1.0e-6                           # [OPT, DEF=1e-10] Aboslute tolerance of the chem. reactor, or pre-factor of the typical values when used
    cvode.solve_type = denseAJ_direct           # [OPT, DEF=GMRES] Linear solver employed for CVODE Newton direction
    cvode.max_order  = 4                        # [OPT, DEF=2] Maximum order of the BDF method in CVODE

Note that the last four parameters belong to the Reactor class of PelePhysics but are specified here for completeness. In particular, CVODE is the adequate choice of integrator to tackle PeleLMeX large time step sizes. Several linear solvers are available depending on whether or not GPU are employed: on CPU, `dense_direct` is a finite-difference direct solver, `denseAJ_direct` is an analytical-jacobian direct solver (preferred choice), `sparse_direct` is an analytical-jacobian sparse direct solver based on the KLU library and `GMRES` is a matrix-free iterative solver; on GPU `GMRES` is a matrix-free iterative solver (available on all the platforms), `sparse_direct` is a batched block-sparse direct solve based on NVIDIA's cuSparse (only with CUDA), `magma_direct` is a batched block-dense direct solve based on the MAGMA library (available with CUDA and HIP.


Linear solvers
--------------

Linear solvers are a key component of PeleLMeX algorithm, separate controls are dedicated to the various solver (MAC projection, nodal projection, diffusion, ...)

::

    #-------------------------LINEAR SOLVERS-----------------------
    nodal_proj.verbose = 1                      # [OPT, DEF=0] Verbose of the nodal projector
    nodal_proj.rtol = 1.0e-11                   # [OPT, DEF=1e-11] Relative tolerance of the nodal projection
    nodal_proj.atol = 1.0e-12                   # [OPT, DEF=1e-14] Absolute tolerance of the nodal projection
    nodal_proj.mg_max_coarsening_level = 5      # [OPT, DEF=100] Maximum number of MG levels (useful when using EB)
    
    mac_proj.verbose = 1                        # [OPT, DEF=0] Verbose of the MAC projector
    mac_proj.rtol = 1.0e-11                     # [OPT, DEF=1e-11] Relative tolerance of the MAC projection
    mac_proj.atol = 1.0e-12                     # [OPT, DEF=1e-14] Absolute tolerance of the MAC projection
    mac_proj.mg_max_coarsening_level = 5        # [OPT, DEF=100] Maximum number of MG levels (useful when using EB)

    diffusion.verbose = 1                       # [OPT, DEF=0] Verbose of the scalar diffusion solve
    diffusion.rtol = 1.0e-11                    # [OPT, DEF=1e-11] Relative tolerance of the scalar diffusion solve
    diffusion.atol = 1.0e-12                    # [OPT, DEF=1e-14] Absolute tolerance of the scalar diffusion solve

    tensor_diffusion.verbose = 1                # [OPT, DEF=0] Verbose of the velocity tensor diffusion solve
    tensor_diffusion.rtol = 1.0e-11             # [OPT, DEF=1e-11] Relative tolerance of the velocity tensor diffusion solve
    tensor_diffusion.atol = 1.0e-12             # [OPT, DEF=1e-14] Absolute tolerance of the velocity tensor diffusion solve

Run-time diagnostics
--------------------

PeleLMeX provides a few diagnostics to check you simulations while it is running as well as adding basic analysis ingredients.

It is often usefull to have an estimate of integrated quantities (kinetic energy, heat release rate, ,..), state extremas
or other overall balance information to get a sense of the status and sanity of the simulation. To this end, it is possible
to activate `temporal` diagnostics performing these reductions at given intervals:

::

    #-------------------------TEMPORALS---------------------------
    peleLM.do_temporals = 1                     # [OPT, DEF=0] Activate temporal diagnostics
    peleLM.temporal_int = 10                    # [OPT, DEF=5] Temporal freq.
    peleLM.do_extremas = 1                      # [OPT, DEF=0] Trigger extremas, if temporals activated
    peleLM.do_mass_balance = 1                  # [OPT, DEF=0] Compute mass balance, if temporals activated
    peleLM.do_species_balance = 1               # [OPT, DEF=0] Compute species mass balance, if temporals activated

The `do_temporal` flag will trigger the creation of a `temporals` folder in your run directory and the following entries 
will be appended to an ASCII `temporals/tempState` file: step, time, dt, kin. energy integral, enstrophy integral, mean pressure
, fuel consumption rate integral, heat release rate integral. Additionnally, if the `do_temporal` flag is activated, one can
turn on state extremas (stored in `temporals/tempExtremas` as min/max for each state entry), mass balance (stored in
`temporals/tempMass`) computing the total mass, dMdt and advective mass fluxes across the domain boundaries as well as the error in
the balance (dMdt - sum of fluxes), and species balance (stored in `temporals/tempSpec`) computing each species total mass, dM_Ydt,
advective \& diffusive fluxes across the domain boundaries, consumption rate integral and the error (dMdt - sum of fluxes - reaction).

Combustion diagnostics often involve the use of a mixture fraction and/or a progress variable, both of which can be defined
at run time and added to the derived variables included in the plotfile. If `mixture_fraction` or `progress_variable` is
added to the `amr.derive_plot_vars` list, one need to provide input for defining those. The mixture fraction is based on
Bilger's element definition and one needs to provide the composition of the 'fuel' and 'oxidizer' tanks using a Cantera-like
format (<species>:<value>) which assumes unspecified species at zero, or a list of floats, in which case all the species must
be specified in the order they appear in the mechanism file.
The progress variable definition in based on a linear combination of the species mass fractions and temperature, and can be
specified in a manner similar to the mixture fraction, providing a list of weights and the prescription of a 'cold' and 'hot'
state:

::

    # ------------------- INPUTS DERIVED DIAGS ------------------
    peleLM.fuel_name = CH4 
    peleLM.mixtureFraction.format = Cantera
    peleLM.mixtureFraction.type   = mass
    peleLM.mixtureFraction.oxidTank = O2:0.233 N2:0.767
    peleLM.mixtureFraction.fuelTank = H2:0.5 CH4:0.5
    peleLM.progressVariable.format = Cantera
    peleLM.progressVariable.weights = CO:1.0 CO2:1.0
    peleLM.progressVariable.coldState = CO:0.0 CO2:0.0
    peleLM.progressVariable.hotState = CO:0.000002 CO2:0.0666


A set of diagnostics available at runtime are currently under development. The following provide an example for extracting
the state variables on a 'x','y' or 'z' aligned plane and writting a 2D plotfile compatible with Amrvis, Paraview or yt:

::

    #--------------------------DIAGNOSTICS------------------------
    
    peleLM.diagnostics = xnormal ynormal
    peleLM.xnormal.type = DiagFramePlane
    peleLM.xnormal.file = xNorm5mm
    peleLM.xnormal.normal = 0
    peleLM.xnormal.center = 0.005
    peleLM.xnormal.int    = 5
    peleLM.xnormal.interpolation = Linear
    
    peleLM.ynormal.type = DiagFramePlane
    peleLM.ynormal.file = yNormCent
    peleLM.ynormal.normal = 1
    peleLM.ynormal.center = 0.0
    peleLM.ynormal.int    = 10
    peleLM.ynormal.interpolation = Quadratic


Run-time control
--------------------

Following some of AMReX's AmrLevel class implementation, PeleLMeX provides a couple of triggers to interact with the code while
it is running. This can be done by adding an empty file to the folder where the simulation is currently running using for 
example:

::

    touch plt_and_continue

The list of available triggers is:

.. list-table:: PeleLMeX run-time triggers
    :widths: 50 100
    :header-rows: 1

    * - File
      - Function
    * - plt_and_continue
      - Write a pltfile to disk and pursue the simulation
    * - chk_and_continue
      - Write a chkfile to disk and pursue the simulation
    * - dump_and_stop
      - Write both pltfile and chkfile to disk and stop the simulation

By default, the code checks if these files exist every 10 time steps, but the user can either increase or decrease the 
frequency using:

::

    amr.message_int      = 20                # [OPT, DEF=10] Frequency for checking the presence of trigger files
    
