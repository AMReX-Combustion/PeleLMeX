Troubleshooting
===============

This section is intended to gather information on common failure mode of PeleLMeX.
Additional information can be found in GitHub issues of `PeleLM <https://github.com/AMReX-Combustion/PeleLM/issues>`_ and `PeleLMeX <https://github.com/AMReX-Combustion/PeleLMeX/issues>`_


Linear solver failure
---------------------

The PeleLMeX algorithm involves multiple linear solves to handle projections
and implicit diffusion. In the event of the solver is enable to solve the
problem, the code will abort with the following message:

```
amrex::Abort::0::MLMG failed !!!
```

or

```
amrex::Abort::0::MLMG failing so lets stop here !!!
```

appearing multiple times when using more than one MPI rank. The first thing to do
is to identify which linear solve is failing and how. To do so, one needs to increase
PeleLMeX, as well as the projection solves verbose (see the `Control <https://amrex-combustion.github.io/PeleLMeX/LMeXControls.html>`_
section for more details on LMeX controls):

::

    peleLM.verbose = 3
    nodal_proj.verbose = 2
    mac_proj.verbose = 2

Note that we focused on the projection solves here because they are generally more
prone to failure than the diffusion ones. You can then restart the simulation
again and identify if the code is failing in the nodal projection, either during the
initial projection (following *Initial velocity projection*) or during the time step
one (following *- oneSDC()::ScalarReaction()  -->*), or in the MAC-projection (right after
*SDC iter [1]*). Then, the linear solver verbose is useful to understand how the solver
fails. If the solver hangs around a small value following an initial reduction of the
residual:

::

    MLMG: # of AMR levels: 1
      # of MG levels on the coarsest AMR level: 9
    MLMG: Initial rhs               = 2666.243975
    MLMG: Initial residual (resid0) = 2666.243975
    MLMG: Iteration   1 Fine resid/bnorm = 0.03858916872
    MLMG: Iteration   2 Fine resid/bnorm = 0.001142880258
    MLMG: Iteration   3 Fine resid/bnorm = 3.300053779e-04
    MLMG: Iteration   4 Fine resid/bnorm = 9.433906375e-06
    MLMG: Iteration   5 Fine resid/bnorm = 2.665697369e-07
    MLMG: Iteration   6 Fine resid/bnorm = 7.40910596e-09
    MLMG: Iteration   7 Fine resid/bnorm = 2.071981144e-10
    MLMG: Iteration   8 Fine resid/bnorm = 2.66772528e-11
    MLMG: Iteration   9 Fine resid/bnorm = 2.568558082e-11
    MLMG: Iteration  10 Fine resid/bnorm = 2.713587827e-11
    MLMG: Iteration  11 Fine resid/bnorm = 2.490776046e-11
    MLMG: Iteration  12 Fine resid/bnorm = 2.41198728e-11
    MLMG: Iteration  13 Fine resid/bnorm = 2.527429436e-11
    MLMG: Iteration  14 Fine resid/bnorm = 2.431036667e-11
    MLMG: Iteration  15 Fine resid/bnorm = 2.479456555e-11
    MLMG: Iteration  16 Fine resid/bnorm = 2.28960372e-11
    MLMG: Iteration  17 Fine resid/bnorm = 2.541484652e-11
    MLMG: Iteration  18 Fine resid/bnorm = 2.522691579e-11
    MLMG: Iteration  19 Fine resid/bnorm = 2.508988366e-11
    ...

it generally means that the required solver tolerance is too small for the problem. The
default relative tolerances of all solvers in PeleLMeX is `1e-11`, but increasing the
resolution, using a small `amr.blocking_factor` (<16) or large flow divergence across
coarse-fine interfaces can lead to the example above. In this case, one can increase the
tolerance of the faulty solver using one of:

::

    nodal_proj.rtol = 5e-11
    mac_proj.rtol   = 5e-11
    diffusion.rtol  = 5e-11

It is sometimes necessary to increase the tolerance up `5e-10`. If you need to go higher
than this ballpark value, it probably indicates that something is wrong in the problem
setup and one should take a closer look at the solution to understand the problem.
Alternatively, the solver can fail as follows:

::

    MLMG: # of AMR levels: 2
      # of MG levels on the coarsest AMR level: 6
    MLMG: Initial rhs               = 395786.0963
    MLMG: Initial residual (resid0) = 395786.0963
    MLMG: Iteration   1 Fine resid/bnorm = 0.009458721163
    MLMG: Iteration   2 Fine resid/bnorm = 1046166408
    MLMG: Iteration   3 Fine resid/bnorm = 5.420966957e+23


In this case, the solver diverges and it is generally a clear indication that the problem
is not properly setup.


Chemistry integration failure
-----------------------------

PeleLMeX relies on `Sundials CVODE <https://computing.llnl.gov/projects/sundials/cvode>`_ to
integrate the stiff ODE resulting of the chemical system (along with advection/diffusion
forcing). CVODE has multiple failure modes, but the most common one appearing in PeleLMeX
will promp a message similar to one of the following:

::

    From CVODE: At t = 0 and h = 6.01889e-195, the corrector convergence test failed repeatedly or with |h| = hmin.```
    From CVODE: At t = 2.459e-6 and h = 6.01889e-16, the corrector convergence test failed repeatedly or with |h| = hmin.```
    [CVODE ERROR]  CVode
        At t = 5.09606e-09, mxstep steps taken before reaching tout.

All of which indicate that the internal sub-stepping algorithm of CVODE did not managed to integrate
the system of ODEs up to the CFL-constrained time step requested by PeleLMeX because CVODE logic
reauired awfully small substep size.

In the case of the first message, one can see that CVODE failed right away (`At t = 0`) which suggests
that the state given to CVODE was wrong. If this happens right at the start of the simulation, your
initial solution is most likely erroneous.

In the case of the second message, the system was integrated up to 2.459e-6 s, but CVODE was not able
to proceed any further as its internal step size dropped to a small value. This could indicates that your
CFL condition is too loose and the chemical stifness can't be properly handled by
CVODE. You can consider reduce your CFL number:

::

    peleLM.cfl = 0.1

if your CFL step size is too large (generally >1e-5 s). e.g. as for a slow, laminar case. This message
can also appear if your state contains species mass fraction undershoots due to poor spatial resolution.
In this case, one can use the following option:

::

    ode.clean_init_massfrac = 1

where the ODE integration is then computed as an increment where the initial species mass fractions
[0-1] bounds are enforced.

