PeleLMeX controls
------------------------

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.
This file needs to specified along with the executable as an `argv` option, for example:

::
    mpirun -np 64 ./PeleLMeX2d.xxx.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inputs file, e.g.:

::
    mpirun -np 64 ./PeleLMeX2d.xxx.ex inputs amr.max_level=2

The available options are divided into groups: those that control primarily AMReX are prefaced with `amr.`, those that are specific to the PeleLM are prefaced by `peleLM.`, while those corresponding to the various pieces of the algorithm are prefaced with specific keys, such that `diffusion`, `nodal_proj`, ... as described below.

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

    #-------------------------AMR CONTROL--------------------------
    amr.n_cell          = 64 64 128        # Level 0 number of cells in each direction   
    amr.v               = 1                # AMR verbose
    amr.max_level       = 1                # maximum level number allowed
    amr.ref_ratio       = 2 2 2 2          # refinement ratio, one for each refinement level
    amr.regrid_int      = 5                # how often to regrid
    amr.n_error_buf     = 1 1 2 2          # number of buffer cells in error est
    amr.grid_eff        = 0.7              # what constitutes an efficient grid
    amr.blocking_factor = 16               # block factor in grid generation (min box size)
    amr.max_grid_size   = 64               # max box size
