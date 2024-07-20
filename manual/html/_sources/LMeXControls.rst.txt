PeleLMeX controls
=================

.. _sec:control:

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.
This file needs to specified along with the executable as an `argv` option, for example:

::

    mpirun -np 64 ./PeleLMeX2d.xxx.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inputs file, e.g.:

::

    mpirun -np 64 ./PeleLMeX2d.xxx.ex inputs amr.max_level=2

The available options are divided into groups: those that control primarily AMReX are prefaced with `amr.`, those that are specific to the PeleLMeX are prefaced by `peleLM.`, while those corresponding to the various pieces of the algorithm are prefaced with specific keys, such that `diffusion`, `nodal_proj`, ... as described below.

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

If specifying boundaries as ``Inflow``, the bcnormal function must be defined
in the ``pelelmex_prob.H`` file for the case to define the inflow conditions.
``Inflow`` boundaries may also be augmented with spatially and temporally
varying turbulent fluctuations using the ``TurbInflow`` utility from
PelePhysics. See the ``Exec/RegTests/TurbInflow`` test for an example of how
to use this capability and the
`PelePhysics  documentation <https://amrex-combustion.github.io/PelePhysics/Utility.html#turbulent-inflows>`_
for the relevant input file flags.

Grid/AMR parameters
-------------------

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

    peleLM.max_grid_size_chem = 32         # [OPT, DEF="None"] Max box size for the Chemistry BoxArray

Load balancing
--------------

*PeleLMeX* relies on two distribution maps (see :doc:`Implementation` for more details).

::

    peleLM.do_load_balancing = 1                    # [OPT, DEF=0] Activate load balancing
    peleLM.load_balancing_method = sfc              # [OPT, DEF="sfc"] AmrCore dmap load balancing method
    peleLM.load_balancing_cost_estimate = ncell     # [OPT, DEF="ncell"] AmrCore dmap balancing cost
    peleLM.chem_load_balancing_method = knapsack    # [OPT, DEF="knapsack"] Chemistry dmap load balancing method
    peleLM.chem_load_balancing_cost_estimate = chemfunctcall_sum # [OPT, DEF="chemfunctcall_sum"] Chemistry dmap balancing cost
    peleLM.load_balancing_efficiency_threshold = 1.05  # What constitute a better dmap ?

The balancing method can be one of `sfc`, `roundrobin` or `knapsack`, while the cost estimate can be one of
`ncell`, `chemfunctcall_avg`, `chemfunctcall_max`, `chemfunctcall_sum`, `userdefined_avg` or `userdefined_sum`. When
using either of the last to option, the user must provide a definition for the `derUserDefined`. If multiple components
are defined in the `derUserDefined` function, the first one is used for load balancing.

Time stepping parameters
------------------------

::

    #-------------------------TIME STEPPING------------------------
    amr.max_step      = 20                 # Maximum number of steps
    amr.stop_time     = 0.001              # Maximum simulation time [s]
    amr.max_wall_time = 1.0                # Maximum wall clock time [hr]
    amr.cfl           = 0.5                # [OPT, DEF=0.7] CFL for advection-controlled dt estimate
    amr.fixed_dt      = 1e-6               # [OPT] optional fixed dt (override CFL condition)
    amr.min_dt        = 1e-11              # [OPT, DEF=1e-12] small time step size limit triggering simulation termination
    amr.init_dt       = 1e-6               # [OPT] optional initial dt (override CFL condition upon initialization)
    amr.dt_shrink     = 0.0001             # [OPT, DEF=1.0] dt factor upon initialization
    amr.dt_change_max = 1.1                # [OPT, DEF=1.1] maximum dt change between consecutive steps

.. note::
   Note that one of `amr.max_step`, `amr.stop_time`, or `amr.max_wall_time` is required, and if more than one is specified,
   the first stopping criterion encountered will lead to termination of the simulation.

IO parameters
-------------

::

    #--------------------------IO CONTROL--------------------------
    amr.plot_int         = 20              # [OPT, DEF=-1] Frequency (as step #) for writing plot file
    amr.plot_per         = 0.002           # [OPT, DEF=-1] Period (time in s) for writing plot file
    amr.plot_per_exact   = 1               # [OPT, DEF=0] Flag to enforce exactly plt_per by shortening dt
    amr.plot_file        = "plt_"          # [OPT, DEF="plt_"] Plot file prefix
    amr.check_int        = 100             # [OPT, DEF=-1] Frequency (as step #) for writing checkpoint file
    amr.check_per        = 0.05            # [OPT, DEF=-1] Period (time in s) for writing checkpoint file
    amr.check_file       = "chk"           # [OPT, DEF="chk"] Checkpoint file prefix
    amr.file_stepDigits  = 6               # [OPT, DEF=5] Number of digits when adding nsteps to plt and chk names
    amr.derive_plot_vars = avg_pressure ...# [OPT, DEF=""] List of derived variable included in the plot files
    amr.plot_speciesState = 0              # [OPT, DEF=0] Force adding state rhoYs to the plot files

    amr.restart          = chk00100        # [OPT, DEF=""] Checkpoint from which to restart the simulation
    amr.initDataPlt      = plt01000        # [OPT, DEF=""] Provide a plotfile from which to extract initial data
    amr.regrid_on_restart = 1              # [OPT, DEF="0"] Trigger a regrid after the data from checkpoint are loaded
    amr.n_files          = 64              # [OPT, DEF="min(256,NProcs)"] Number of files to write per level
    
Refinement controls
-------------------

Refinement in PeleLMeX is controlled by a set of 'Tagging' criterion listed under the `amr.refinement_indicators`
key. For each criteriq, the user needs to supply a definition. For example, the following provides a complete
overview of the available controls:

::

    amr.refinement_indicators gthan lthan adjd box1

    amr.gthan.max_level     = 3
    amr.gthan.value_greater = 0.005
    amr.gthan.field_name    = x_velocity

    amr.lthan.max_level     = 4
    amr.lthan.value_less    = 400.0
    amr.lthan.field_name    = temp
    amr.lthan.start_time    = 0.001
    amr.lthan.end_time      = 0.005

    amr.adjd.max_level                   = 2
    amr.adjd.adjacent_difference_greater = 0.05
    amr.adjd.field_name                  = density

    amr.box1.max_level      = 1
    amr.box1.in_box_lo      = 0.0 0.0 0.0
    amr.box1.in_box_hi      = 0.01 0.01 0.05

The `field_name` can be any of the state or derived variables (see below) component. Additional controls specific
to embedded boundaries are discussed below.

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
    * - `diffcoeff`
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
      - Vorticity magnitude
    * - `vorticity`
      - AMREX_SPACEDIM*2-3
      - VortZ (2D) or VortX, VortY, VortZ (3D)
    * - `Qcrit`
      - 1
      - Q-Criterion : :math:`0.5(|\boldsymbol{\Omega}|^2 - |\boldsymbol{S}|^2)`
    * - `kinetic_energy`
      - 1
      - Kinetic energy: 0.5 * rho * (u^2+v^2+w^2)
    * - `enstrophy`
      - 1
      - enstrophy: 0.5 * rho * (\omega_x^2+\omega_y^2+\omega_z^2)
    * - `HeatRelease`
      - 1
      - Heat release rate from chem. reactions
    * - `rhominsumrhoY`
      - 1
      - Rho minus sum of rhoYs, for debug purposes
    * - `coordinates`
      - AMREX_SPACEDIM
      - Cell-center coordinates
    * - `DistributionMap`
      - 1
      - The MPI-rank of each box
    * - `derUserDefined`
      - ?
      - A user-defined derived which number of components is provided by the user (see below).

Note that `mixture_fraction` and `progress_variable` requires additional inputs from the users as described below.
The `derUserDefined` allow the user to define its own derived variable which can comprise several components. To do
so, the user need to copy the Source/DeriveUserDefined.cpp file into his run folder and update the file. The number of
components is defined based on the size of the vector returned by pelelmex_setuserderives().

PeleLMeX algorithm
------------------

::

    #-----------------------PELE CONTROL-----------------------
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
    peleLM.do_periodic_channel = 0         # [OPT, DEF= 0] Add an automatic pressure gradient to maintain initial condition mass flow rate in periodic channel
    peleLM.periodic_channel_dir = 2        # [OPT, DEF= -1] Required if do_periodic_channel != 0. Direction to apply pressure gradient.
    peleLM.closed_chamber = 0              # [OPT] Override the automatic detection of closed chamber (based on Outflow(s))
    peleLM.floor_species = 0               # [OPT, DEF=0] Crudely enforce mass fraction positivity
    peleLM.deltaT_verbose = 0              # [OPT, DEF=0] Verbose of the deltaT iterative solve algorithm
    peleLM.deltaT_iterMax = 5              # [OPT, DEF=10] Maximum number of deltaT iterations
    peleLM.deltaT_tol = 1e-10              # [OPT, DEF=1.e-10] Tolerance of the deltaT solve
    peleLM.evaluate_vars =...              # [OPT, DEF=""] In evaluate mode, list unitTest: diffTerm, divU, instRR, transportCC
    peleLM.do_patch_flow_variables = false # [OPT, DEF=false] Enable user-defined flow variable patching after reading a plot solution file

Transport coefficients and LES
------------------------------

::

    #-----------------------DIFFUSION AND LES MODEL CONTROL-----------------------
    peleLM.fixed_Le = 0                    # [OPT, DEF=0] Use a fixed Lewis number approximation for species diffusivities
    peleLM.fixed_Pr = 0                    # [OPT, DEF=0] Use a fixed Prandtl number approximation for thermal diffusivity
    peleLM.Prandtl = 0.7                   # [OPT, DEF=0.7] If fixed_Pr or doing LES, specifies the Prandtl number
    peleLM.Schmidt = 0.7                   # [OPT, DEF=0.7] If doing LES, specifies the Schmidt number
    peleLM.Lewis = 1.0                     # [OPT, DEF=1.0] If fixed_Le, specifies the Lewis number

    peleLM.les_model = "None"              # [OPT, DEF="None"] Model to compute turbulent viscosity: None, Smagorinsky, WALE, Sigma
    peleLM.les_cs_smag = 0.18              # [OPT, DEF=0.18] If using Smagorinsky LES model, provides model coefficient
    peleLM.les_cm_wale = 0.60              # [OPT, DEF=0.60] If using WALE LES model, provides model coefficient
    peleLM.les_cs_sigma = 1.35             # [OPT, DEF=1.35] If using Sigma LES model, provides model coefficient
    peleLM.les_v = 0                       # [OPT, DEF=0] Verbosity level for LES model
    peleLM.plot_les = 0                    # [OPT, DEF=0] If doing LES, whether to plot the turbulent viscosity

Chemistry integrator
--------------------

::

    #-----------------------CHEMISTRY CONTROL----------------------
    peleLM.chem_integrator   = "ReactorCvode"   # Chemistry integrator, from PelePhysics available list
    peleLM.use_typ_vals_chem = 1                # [OPT, DEF=1] Use Typical values to scale components in the reactors
    peleLM.typical_values_reset_int = 5         # [OPT, DEF=10] Frequency at which the typical values are updated
    ode.rtol = 1.0e-6                           # [OPT, DEF=1e-10] Relative tolerance of the chem. reactor
    ode.atol = 1.0e-6                           # [OPT, DEF=1e-10] Absolute tolerance of the chem. reactor, or pre-factor of the typical values when used
    cvode.solve_type = denseAJ_direct           # [OPT, DEF=GMRES] Linear solver employed for CVODE Newton direction
    cvode.max_order  = 4                        # [OPT, DEF=2] Maximum order of the BDF method in CVODE
    cvode.max_substeps = 10000                  # [OPT, DEF=10000] Maximum number of substeps for the linear solver in CVODE

Note that the last five parameters belong to the Reactor class of PelePhysics but are specified here for completeness. In particular, CVODE is the adequate choice of integrator to tackle PeleLMeX large time step sizes. Several linear solvers are available depending on whether or not GPU are employed: on CPU, `dense_direct` is a finite-difference direct solver, `denseAJ_direct` is an analytical-jacobian direct solver (preferred choice), `sparse_direct` is an analytical-jacobian sparse direct solver based on the KLU library and `GMRES` is a matrix-free iterative solver; on GPU `GMRES` is a matrix-free iterative solver (available on all the platforms), `sparse_direct` is a batched block-sparse direct solve based on NVIDIA's cuSparse (only with CUDA), `magma_direct` is a batched block-dense direct solve based on the MAGMA library (available with CUDA and HIP. Different `cvode.solve_type` should be tried before increasing the `cvode.max_substeps`.

.. note::
   The default chemistry integrator is 'ReactorNull' which do not include the chemical source terms.

Embedded Geometry
-----------------

`PeleLMeX` geometry relies on AMReX implementation of the EB method. Simple geometrical objects
can thus be constructed using `AMReX internal parser <https://amrex-codes.github.io/amrex/docs_html/EB.html>`_.
For instance, setting up a sphere of radius 5 mm can be achieved:

::

    eb2.geom_type = sphere
    eb2.sphere_radius = 0.005
    eb2.sphere_center = 0.0 0.0 0.0
    eb2.sphere_has_fluid_inside = 0
    eb2.small_volfrac = 1.0e-4
    eb2.maxiter = 200

The `eb2.small_volfrac` controls volume fraction that are deemed too small and eliminated from the EB representation.
This operation is done iteratively and the maximum number of iteration is prescribed by `eb2.maxiter`.
For most applications, a single AMReX object is insufficient to represent the geometry. AMReX enable to combine
objects using constructive solid geometry (CSG) in order to create complex geometry. It is up to users to define
the combination of basic elements leading to their desired geometry. To switch to a user-defined EB definition, one
must set:

::

    eb2.geom_type = UserDefined

and then implement the actual geometry definition in a `EBUserDefined.H` file located in the run folder (and add
to the GNUmakefile using `CEXE_headers += EBUserDefined.H`). An example of such implementation is available in the
``Exec/Case/ChallengeProblem`` folder. Example of more generic EB problems are also found in the ``Exec/RegTest/EB_*``
folders.

In addition to the input keys presented above, a set of `PeleLMeX`-specific keys are available in order to control refinement at the EB:

::

    peleLM.refine_EB_type = Static
    peleLM.refine_EB_max_level = 1
    peleLM.refine_EB_buffer = 2.0

By default, the EB is refined to the `amr.max_level`, which can lead to undesirably high number of cells
close to the EB when the physics of interest might be elsewhere. The above lines enable to limit the
EB-level to level 1 (must be below `amr.max_level`) and a derefinement strategy is adopted to ensure
that fine-grid patches do not cross the EB boundary. The last parameter set a safety margin to increase
how far the derefinement is applied in order to account for grid-patches diagonals and proper nesting constrains.
Note that the parameter do not ensure explicitly coarse-fine/EB crossings are avoided and the code will fail when this happens.

AMReX generates the EB at the finest level (specified by `amr.max_level`) and subsequently coarsen the resulting
EB data to coarser AMR levels (as well as multigrid levels) in order to ensure consistency across levels. As a consequence,
increasing `amr.max_level` during the course of a simulation can lead in small changes to the EB, potentially uncovering
previously covered regions for which the solver cannot provide initial conditions upon restart. It is thus advised to
set `amr.max_level` to the desired maximum level you ever plan to use in the simulation, and limit the level at which
the EB is refined (with `peleLM.refine_EB_max_level`) and refinement criteria are applied (with `amr.*****.max_level`)
to a lower value while starting the simulation. A side effect of this process is that the generation of the EB can
be excessively long if `amr.max_level` is large. To speed up the EB generation, one can use the following keys:

::

    eb2.num_coarsen_opt = 3
    eb2.max_grid_size = 64

Setting `eb2.num_coarsen_opt = 3` effectively speed the EB generation by coarsening the level at which the recursive EB intersection
algorithm starts by a factor 2^3. A large value of this parameter can lead to an erroneous EB representation, missing small features,
and should thus be kept at or below 4. The EB intersection algorithm relies on chopping the domain into boxes and finding those
fully or partially covered. The size of these boxes can be controlled with `eb2.max_grid_size`, and can be adjusted to better
match the number of MPI ranks used in the simulation.


It is also possible to change the default adiabatic EB wall condition to an isothermal EB. To do so, one need to switch the following
flag:

::

    peleLM.isothermal_EB = 1

The user is now responsible for providing the wall temperature *on all the EB walls*, but adiabtic wall can still be specified.
Control over the local EB thermal boundary condition is provided through the `setEBState` and `setEBType` functions, also
defined in the `EBUserDefined.H` already used above to provide a user-defined EB geometry. Example of isothermal EBs are provided
in ``Exec/RegTest/EB_BackwardStepFlame`` and ``Exec/RegTest/EB_FlowPastCylinder`` tests.

.. note::
   Note that when using isothermal EB in combination with LES, the thermal diffusion coefficient employed to compute the EB boundary thermal flux only uses the molecular contribution.

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

Active control
--------------

`PeleLMeX` includes an active control mechanism to enable statistically steady simulations of flames
maintaining the flame at a fixed position in the domain. An example of this feature is provided in
the triple flame tutorial :doc:`Tutorials_TripleFlame`.

.. note::
   To enable active control, a ``FlowControllerData FCData`` object must be added to the problem ``ProbParm``!

During the course of the simulation, the FlowControllerData is updated based on the flame position to allow
the user to set the inflow velocity. The following options are available when using active control: ::

    #---------------------- AC CONTROL -------------------------------
    active_control.on = 1                     # [OPT, DEF=0] Use AC ?
    active_control.use_temp = 1               # [OPT, DEF=1] Default in fuel mass, rather use iso-T position ?
    active_control.temperature = 1400.0       # [OPT, DEF=-1] Value of iso-T ?
    active_control.tau = 5.0e-4               # [OPT, DEF=0.0] Control tau (should ~ 10 dt)
    active_control.height = 0.01              # [OPT, DEF=0.0] Where is the flame held ?
    active_control.v = 1                      # [OPT, DEF=0] verbose
    active_control.method = 1                 # [OPT, DEF=2] Controller: 1 - Linear, 2 - Quadratic, 3 - Weighted quadratic
    active_control.velMax = 2.0               # [OPT, DEF=-1.0] limit inlet velocity, only used when positive
    active_control.changeMax = 0.1            # [OPT, DEF=1.0] limit inlet velocity changes (absolute m/s)
    active_control.flow_dir  = 1              # [OPT, DEF=AMREX_SPACEDIM-1] flame main direction
    active_control.AC_history  = AChist       # [OPT, DEF=AC_history] Control history file, read upon restart
    active_control.npoints_average = 5        # [OPT, DEF=3] Number of previous steps using to estimate new velocity
    active_control.pseudo_gravity = 1         # [OPT, DEF=0] add density proportional force to compensate for the acceleration
                                              #           of the gas due to inlet velocity changes

Run-time diagnostics
--------------------

`PeleLMeX` provides a few diagnostics to check you simulations while it is running as well as adding basic analysis ingredients.

It is often useful to have an estimate of integrated quantities (kinetic energy, heat release rate, ,..), state extremas
or other overall balance information to get a sense of the status and sanity of the simulation. To this end, it is possible
to activate `temporal` diagnostics performing these reductions at given intervals:

::

    #-------------------------TEMPORALS---------------------------
    peleLM.do_temporals = 1                     # [OPT, DEF=0] Activate temporal diagnostics
    peleLM.temporal_int = 10                    # [OPT, DEF=5] Temporal freq.
    peleLM.do_extremas = 1                      # [OPT, DEF=0] Trigger extremas, if temporals activated
    peleLM.do_mass_balance = 1                  # [OPT, DEF=0] Compute mass balance, if temporals activated
    peleLM.do_species_balance = 1               # [OPT, DEF=0] Compute species mass balance, if temporals activated
    peleLM.do_patch_mfr=1                       # [OPT, DEF=0] Activate patch based species flux diagbostics
    peleLM.bpatch.patchnames= <patch_name1 patch_name2 ..> # List of patchnames
    
    bpatch.patch_name1.patchtype=full-boundary             # patchtype one of "full-boundary", "circle, "rectangle", "circle-annular" or "rectangle-annular"	 
    bpatch.patch_name1.boundary_direction=2                # patch normal direction 
    bpatch.patch_name1.boundary_lo_or_hi=0                 # patch in low or high side of boundary    
    bpatch.patch_name1.species= O2 N2                      # list of species names
    
    bpatch.patch_name2.patchtype=circle                    # patchtype one of "full-boundary", "circle, "rectangle", "circle-annular" or "rectangle-annular"	 
    bpatch.patch_name2.boundary_direction=2                # patch normal direction 
    bpatch.patch_name2.boundary_lo_or_hi=0                 # patch in low or high side of boundary   
    bpatch.patch_name2.patch_circle_radius=0.1             # radius of the patch
    bpatch.patch_name2.patch_circle_center=0.0 0.0 0.0     # coordinates of patch center 
    bpatch.patch_name2.species= O2 N2                      # list of species names
    
    bpatch.patch_name3.patchtype=rectangle      	 
    bpatch.patch_name3.boundary_direction=2                # patch normal direction 
    bpatch.patch_name3.boundary_lo_or_hi=0                 # patch in low or high side of boundary   
    bpatch.patch_name3.patch_rectangle_lo=0.0 0.0 0.0      # coordinates of low corner of rectangle
    bpatch.patch_name3.patch_rectangle_hi=1.0 1.0 1.0      # coordinates of high corner of rectangle
    bpatch.patch_name3.species= O2 N2                      # list of species names
    
    bpatch.patch_name4.patchtype=circle-annular  	 
    bpatch.patch_name4.boundary_direction=2                # patch normal direction 
    bpatch.patch_name4.boundary_lo_or_hi=0                 # patch in low or high side of boundary   
    bpatch.patch_name4.patch_circ_ann_center= 0.0 0.0 0.0  # center of annular circle
    bpatch.patch_name4.patch_circ_ann_inner_radius=0.1     # coordinates of patch center
    bpatch.patch_name4.patch_circ_ann_outer_radius=0.2     # coordinates of patch center 
    bpatch.patch_name4.species= O2 N2                      # list of species names
    
    bpatch.patch_name5.patchtype=rectangle-annular         
    bpatch.patch_name5.boundary_direction=2                     # patch normal direction 
    bpatch.patch_name5.boundary_lo_or_hi=0                      # patch in low or high side of boundary   
    bpatch.patch_name5.patch_rect_ann_outer_lo = -1.0 -1.0 -1.0 # coordinates of low corner of outer rectangle 
    bpatch.patch_name5.patch_rect_ann_outer_hi =  1.0  1.0  1.0 # coordinates of high corner of outer rectangle  
    bpatch.patch_name5.patch_rect_ann_inner_lo = -0.5 -0.5 -0.5 # coordinates of low corner of inner rectangle 
    bpatch.patch_name5.patch_rect_ann_inner_hi =  0.5  0.5  0.5 # coordinates of high corner of inner rectangle 
    bpatch.patch_name5.species= O2 N2                           # list of species names

The `do_temporal` flag will trigger the creation of a `temporals` folder in your run directory and the following entries
will be appended to an ASCII `temporals/tempState` file: step, time, dt, kin. energy integral, enstrophy integral, mean pressure
, fuel consumption rate integral, heat release rate integral. Additionally, if the `do_temporal` flag is activated, one can
turn on state extremas (stored in `temporals/tempExtremas` as min/max for each state entry), mass balance (stored in
`temporals/tempMass`) computing the total mass, dMdt and advective mass fluxes across the domain boundaries as well as the error in
the balance (dMdt - sum of fluxes), and species balance (stored in `temporals/tempSpec`) computing each species total mass, dM_Ydt,
advective \& diffusive fluxes across the domain boundaries, consumption rate integral and the error (dMdt - sum of fluxes - reaction).
Users can also monitor species advective fluxes through specific regions of the domain boundaries (called as boundary patches).
Patches can be defined on the low or high sides of non-embedded boundaries through the use of pre-defined shapes such as `circle`,
`rectangle`,`circle-annular`, `rectangle-annular` and `full-boundary`. The zero AMR level, advective fluxes of each of the user-specified species will be 
reported in the ASCII `temppatchmfr` file in the temporals folder. 

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


Analysing the data a-posteriori can become extremely cumbersome when dealing with extreme datasets.
PeleLMeX offers a set of diagnostics available at runtime and more are under development.
Currently, the list of diagnostic contains:

* `DiagFramePlane` : extract a plane aligned in the 'x','y' or 'z' direction across the AMR hierarchy, writing
  a 2D plotfile compatible with Amrvis, Paraview or yt. Only available for 3D simulations.
* `DiagPDF` : extract the PDF of a given variable and write it to an ASCII file.
* `DiagConditional` : extract statistics (average and standard deviation, integral or sum) of a
  set of variables conditioned on the value of given variable and write it to an ASCII file.

When using `DiagPDF` or `DiagConditional`, it is possible to narrow down the diagnostic to a region of interest
by specifying a set of filters, defining a range of interest for a variable. Note also the for these two diagnostics,
fine-covered regions are masked. The following provide examples for each diagnostic:

::

    #--------------------------DIAGNOSTICS------------------------

    peleLM.diagnostics = xnormP condT pdfTest

    peleLM.xnormP.type = DiagFramePlane                             # Diagnostic type
    peleLM.xnormP.file = xNorm5mm                                   # Output file prefix
    peleLM.xnormP.normal = 0                                        # Plane normal (0, 1 or 2 for x, y or z)
    peleLM.xnormP.center = 0.005                                    # Coordinate in the normal direction
    peleLM.xnormP.int    = 5                                        # Frequency (as step #) for performing the diagnostic
    peleLM.xnormP.interpolation = Linear                            # [OPT, DEF=Linear] Interpolation type : Linear or Quadratic
    peleLM.xnormP.field_names = x_velocity mag_vort density         # List of variables outputted to the 2D pltfile
    peleLM.xnormP.n_files = 2                                       # [OPT, DEF="min(256,NProcs)"] Number of files to write per level

    peleLM.condT.type = DiagConditional                             # Diagnostic type
    peleLM.condT.file = condTest                                    # Output file prefix
    peleLM.condT.int  = 5                                           # Frequency (as step #) for performing the diagnostic
    peleLM.condT.filters = xHigh stoich                             # [OPT, DEF=None] List of filters
    peleLM.condT.xHigh.field_name = x                               # Filter field
    peleLM.condT.xHigh.value_greater = 0.006                        # Filter definition : value_greater, value_less, value_inrange
    peleLM.condT.stoich.field_name = mixture_fraction               # Filter field
    peleLM.condT.stoich.value_inrange = 0.053 0.055                 # Filter definition : value_greater, value_less, value_inrange
    peleLM.condT.conditional_type = Average                         # Conditional type : Average, Integral or Sum
    peleLM.condT.nBins = 50                                         # Number of bins for the conditioning variable
    peleLM.condT.condition_field_name = temp                        # Conditioning variable name
    peleLM.condT.field_names = HeatRelease I_R(CH4) I_R(H2)         # List of variables to be treated

    peleLM.pdfTest.type = DiagPDF                                   # Diagnostic type
    peleLM.pdfTest.file = PDFTest                                   # Output file prefix
    peleLM.pdfTest.int  = 5                                         # Frequency (as step #) for performing the diagnostic
    peleLM.pdfTest.filters = innerFlame                             # [OPT, DEF=None] List of filters
    peleLM.pdfTest.innerFlame.field_name = temp                     # Filter field
    peleLM.pdfTest.innerFlame.value_inrange = 450.0 1500.0          # Filter definition : value_greater, value_less, value_inrange
    peleLM.pdfTest.nBins = 50                                       # Number of bins for the PDF
    peleLM.pdfTest.normalized = 1                                   # [OPT, DEF=1] PDF is normalized (i.e. integral is unity) ?
    peleLM.pdfTest.volume_weighted = 1                              # [OPT, DEF=1] Computation of the PDF is volume weighted ?
    peleLM.pdfTest.range = 0.0 2.0                                  # [OPT, DEF=data min/max] Specify the range of the PDF
    peleLM.pdfTest.field_name = x_velocity                          # Variable of interest

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
