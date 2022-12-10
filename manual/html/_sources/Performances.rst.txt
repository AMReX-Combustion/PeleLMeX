Performances
============

`PeleLMeX` development was driven by the need to create a simulation code efficiently
leveraging the computational power of ExaScale super-computers. As mentioned earlier,
`PeleLMeX` is built upon the AMR library `AMReX` and inherits most of its High Performance Computing
features.

`PeleLMeX` parallel paradigm is based on an MPI+`X` appraoch, where `X` can be OpenMP, or any of
CUDA, HIP or SYCL, for Nvidia, AMD and Intel GPUs vendor, respectively. The actual performances
gain of using accelerator within `PeleLMeX` is a moving target as both hardware and software are
continously improving. In the following we demonstrate the gain at a given time (specified and
subject to updates) and on selected platforms.


Single node performances: FlameSheet case
-----------------------------------------

Case description
^^^^^^^^^^^^^^^^

The simple case of a laminar premixed flame with harmonic perturbation can be found in
`Exec/RegTests/FlameSheet`. For the following test, the mixture is composed of
dodecane and air at ambient temperature and pressure. The chemical mechanism used consist
of 35 transported species and 18 species assumed in Quqsi-Steady State (QSS) and the `Simple`
transport model with the `Fuego` EOS is used:

::

    Chemistry_Model = dodecane_lu_qss
    Eos_Model       = Fuego
    Transport_Model = Simple

The initial solution is provided from a Cantera simulation and averaged on the cartesian grid.
The input file `Exec/RegTests/FlameSheet/inputs.3d_DodecaneQSS` is used, with modifications detailed hereafter.
Simulations are conducted at a fixed time step size for 16 steps, bypassing the initial reduction of
the time step size usually employed to remove artifacts from the initial data:

::

    amr.max_step = 16
    amr.dt_shrink = 1.0
    amr.fixed_dt = 2.5e-7

Additionnaly, all the tests on GPUs are conducted using the MAGMA dense-direct solver to solve for
the Newton direction within CVODE's non-linear integration.

::

    cvode.solve_type = magma_direct

and the dense direct analytical Jacobian solver on CPUs:

::

    cvode.solve_type = denseAJ_direct


The actual number of cells in each direction and the number of levels will depends on the amount
of memory available on the different platform and will be specified later on.

Results on Perlmutter (NERSC)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Perlmutter's `GPU nodes <https://docs.nersc.gov/systems/perlmutter/architecture/#gpu-nodes>`_ consists of a single AMD EPYC 7763 (Milan)
CPU connected to 4 NVIDIA A100 GPUs. The `CPU nodes <https://docs.nersc.gov/systems/perlmutter/architecture/#cpu-nodes>`_ consists of
two of the same AMD EPYC, 64-cores CPUs. When running on the GPU node, `PeleLMeX` will use 4 MPI ranks with each access to one A100, while
when running on a CPU node, we will use 128 MPI-ranks.

The FlameSheet case is ran using 2 levels of refinement (3 levels total) and the following domain size and cell count:

::

    geometry.prob_lo     = 0.0 0.0 0.0        # x_lo y_lo (z_lo)
    geometry.prob_hi     = 0.008 0.016 0.016  # x_hi y_hi (z_hi)

    amr.n_cell           = 32 64 64
    amr.max_level        = 2

leading to an initial cell count of 3.276 M, i.e. 0.8M/cells per GPU. The git hashes of `PeleLMeX` and its dependencies for
these tests are:

::

     ================= Build infos =================
     PeleLMeX    git hash: v22.12-dirty
     AMReX       git hash: 22.12-1-g4a53367b1-dirty
     PelePhysics git hash: v0.1-1052-g234a8089-dirty
     AMReX-Hydro git hash: d959ee9
     ===============================================

The graph below compares the timings of the two runs obtained from `AMReX` TinyProfiler.
Inclusive averaged data are presented here, for separates portion of the `PeleLMeX` algorithm
(see the `algorithm page <https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#pelelmex-algorithm>`_ for more
details):


.. figure:: images/performances/PMF/SingleNodePMF_PM.png
   :align: center
   :figwidth: 90%

The total time comparison shows close to a 4x speed-up on a node basis on this platform, with the AMD Milan CPU being amongst
the most performant to date. The detailed distribution of the computational time within each run highlight the dominant contribution
of the stiff chemistry integration, especially on the GPU.
