.. role:: cpp(code)
   :language: c++

.. _sec:insitu:

In-Situ Visualization with Ascent
==================================

`PeleLMeX` supports in-situ visualization via `Ascent
<https://ascent.readthedocs.io>`_, an open-source many-core capable
lightweight in-situ visualization and analysis library developed as part of the
`Alpine <https://alpine-dav.github.io/ascent/>`_ project. Ascent uses `Conduit
<https://llnl-conduit.readthedocs.io>`_ to describe and pass simulation data,
and the `Viskores <https://github.com/Viskores/viskores>`_ library for
rendering on both CPU and GPU. Since the solver state is passed directly to
Ascent without writing to disk, in-situ rendering eliminates the I/O bottleneck
of traditional post-hoc workflows and is well-suited to large-scale GPU runs.

The PeleLMeX Ascent integration publishes exactly the same fields as
``WritePlotFile()``, controlled at runtime by the same input file flags.
Any field visible in a plotfile is also available for in-situ rendering.

.. note::
   Ascent and Conduit must be built and installed before enabling in-situ
   visualization. Refer to the `Ascent build documentation
   <https://ascent.readthedocs.io/en/latest/BuildingAscent.html>`_ and the
   `Conduit build documentation
   <https://llnl-conduit.readthedocs.io/en/latest/building.html>`_ for
   instructions. The ``build_ascent.sh`` script provided in the Ascent
   repository (``scripts/build_ascent/build_ascent.sh``) builds both Ascent
   and Conduit together and is the recommended approach.

.. _sec:insitu::build:

Building with Ascent
--------------------

Ascent support is enabled at compile time by passing the following flags to
the GNUmake build system. AMReX's GNUmake locates MPI via the compiler
wrappers in ``PATH``, so ensure your MPI installation's ``bin/`` directory
is in ``PATH`` before running ``make``: ::

    make -j8 USE_CUDA=TRUE CUDA_ARCH=89         \
             USE_ASCENT=TRUE                    \
             ASCENT_DIR=/path/to/ascent/install \
             USE_CONDUIT=TRUE                   \
             CONDUIT_DIR=/path/to/conduit/install

Replace ``/path/to/ascent/install`` and ``/path/to/conduit/install`` with the
paths to your Ascent and Conduit installations. ``CUDA_ARCH`` should match your
GPU's compute capability (e.g., ``89`` for NVIDIA RTX 40-series, ``80`` for
A100). Omit the ``USE_CUDA`` flags to build a CPU-only Ascent-enabled
executable.

These flags can be combined with any other physics flags (``USE_SOOT``,
``USE_RADIATION``, ``USE_PARTICLES``, ``USE_PLASMA``). The Ascent integration
automatically publishes the additional fields for each compiled physics module
when the corresponding runtime flags are active.

.. warning::
   Ascent, Conduit, and PeleLMeX must all be built against the **same MPI
   installation**. Building any component against a different MPI will cause
   ABI mismatches and runtime failures when ``libascent_mpi.so`` or
   ``libconduit_mpi.so`` are loaded. Verify before building: ::

       which mpicc        # must point to your chosen MPI installation
       mpicc --version    # confirm the version matches across all components

   The MPI implementation itself does not matter — OpenMPI, MPICH, MVAPICH,
   Intel MPI, and Cray MPI are all supported — but the same installation must
   be used throughout. The ``LD_LIBRARY_PATH`` must also include the MPI, Ascent,
   and Conduit library directories at runtime.

.. _sec:insitu::runtime:

Runtime activation
------------------

Ascent is activated at runtime by adding the following to the input file: ::

    ascent.plot_int = 10    # call Ascent every 10 time steps

When ``ascent.plot_int`` is not set or is negative, Ascent is fully disabled
with no runtime overhead.

Ascent reads two YAML files from the run directory automatically:

- ``ascent_actions.yaml`` — defines what to render (scenes, pipelines, filters).
  This file is **required** for Ascent to produce any output.
- ``ascent_options.yaml`` — optional runtime configuration. The most common use
  is to override the rendering backend: ::

    runtime:
      viskores:
        backend: openmp    # valid values: cuda, openmp, serial, kokkos

  When ``ascent_options.yaml`` is absent or no backend is specified, Ascent
  selects the highest-performance backend available in your build, using the
  priority order CUDA → OpenMP → Kokkos → Serial. Which backends are available
  depends on the flags passed to ``build_ascent.sh``: a build with
  ``enable_cuda=ON`` will default to CUDA; a CPU-only build with
  ``enable_openmp=ON`` will default to OpenMP. Override explicitly when you
  want to free GPU memory for the solver (``backend: openmp``) or force a
  specific device in a multi-GPU environment.

For a full reference of available Ascent actions (contours, volume rendering,
Cinema databases, triggers, expressions, and more), see the
`Ascent actions documentation
<https://ascent.readthedocs.io/en/latest/Actions/Actions.html>`_.

.. _sec:insitu::fields:

Published fields
----------------

PeleLMeX passes the same field set to Ascent as it writes to plotfiles.
The tables below list every field available for rendering in
``ascent_actions.yaml``, grouped by the runtime flag or compile-time option
that controls their inclusion. Fields are always referred to by their exact
string name in the yaml ``field:`` key.

Base state
^^^^^^^^^^

Always published. Species fields are controlled by ``amr.plot_speciesState``.

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name(s)
     - Components
     - Description
   * - ``x_velocity``, ``y_velocity`` [, ``z_velocity``]
     - SPACEDIM
     - Velocity components
   * - ``density``
     - 1
     - Mixture density :math:`\rho`
   * - ``rho.Y(<species>)``
     - NUM_SPECIES
     - Species partial densities :math:`\rho Y_k`. Published when
       ``amr.plot_speciesState = 1`` (default). Set to ``0`` to suppress.
   * - ``rhoh``
     - 1
     - Mixture enthalpy :math:`\rho h`
   * - ``temp``
     - 1
     - Temperature
   * - ``RhoRT``
     - 1
     - :math:`\rho R T` (thermodynamic pressure proxy)
   * - ``divu``
     - 1
     - Velocity divergence constraint (when ``peleLM.has_divu = 1``)
   * - ``gradpx`` [, ``gradpy``, ``gradpz``]
     - SPACEDIM
     - Pressure gradient components (when ``peleLM.plot_grad_p = 1``)

Reaction rates
^^^^^^^^^^^^^^

Published when ``peleLM.do_react = 1`` (reacting flow) and
``peleLM.plot_react = 1`` (default when reacting).

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name(s)
     - Components
     - Description
   * - ``I_R(<species>)``
     - NUM_SPECIES
     - Species reaction rates :math:`\dot{\omega}_k`
   * - ``FunctCall``
     - 1
     - CVODE integrator function call count per cell (chemistry stiffness diagnostic)
   * - ``HeatRelease``
     - 1
     - Heat release rate. Published when additionally ``peleLM.plot_heatRelease = 1``.

Derived variables
^^^^^^^^^^^^^^^^^

Published for each name listed in ``amr.derive_plot_vars``. The full list of
available derived variables and their descriptions is given in the
:doc:`LMeXControls` page under *PeleLMeX derived variables*. Commonly useful
for in-situ visualization:

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``avg_pressure``
     - 1
     - Cell-averaged pressure from nodal :math:`\pi`
   * - ``mag_vort``
     - 1
     - Vorticity magnitude :math:`|\boldsymbol{\omega}|`
   * - ``vorticity``
     - 1 (2D) / 3 (3D)
     - Vorticity components
   * - ``Qcrit``
     - 1
     - Q-criterion
   * - ``kinetic_energy``
     - 1
     - :math:`\frac{1}{2} \rho |\mathbf{u}|^2`
   * - ``enstrophy``
     - 1
     - :math:`\frac{1}{2} \rho |\boldsymbol{\omega}|^2`
   * - ``viscosity``
     - 1
     - Mixture dynamic viscosity
   * - ``mass_fractions``
     - NUM_SPECIES
     - Species mass fractions ``Y(<species>)``
   * - ``mixture_fraction``
     - 1
     - Bilger mixture fraction (requires additional inputs, see :doc:`LMeXControls`)
   * - ``progress_variable``
     - 1
     - Progress variable (requires additional inputs, see :doc:`LMeXControls`)

LES turbulent viscosity
^^^^^^^^^^^^^^^^^^^^^^^

Published when ``peleLM.les_model`` is set to a non-``None`` model **and**
``peleLM.plot_les = 1``. Computing the turbulent viscosity at the plot time
requires a full velocity gradient tensor evaluation; set ``peleLM.plot_les = 0``
to suppress this cost while retaining the resolved LES flow fields above.

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``viscturb``
     - 1
     - Turbulent (SGS) viscosity, face-to-cell interpolated

Soot (``PELE_USE_SOOT``)
^^^^^^^^^^^^^^^^^^^^^^^^

Published when compiled with ``USE_SOOT=TRUE`` and ``peleLM.do_soot_solve = 1``.
Field names are assigned by the soot model from the moment indices; the exact
names are mechanism-dependent and can be discovered at runtime (see
:ref:`sec:insitu::discovery`).

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``soot_N``, ``soot_N0``
     - 1 each
     - Total and nucleation soot number densities
   * - ``soot_S``
     - 1
     - Soot surface area density
   * - ``soot_fv``
     - 1
     - Soot volume fraction

Radiation (``PELE_USE_RADIATION``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Published when compiled with ``USE_RADIATION=TRUE`` and
``peleLM.do_rad_solve = 1``.

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``rad.G``
     - 1
     - Mean radiative intensity :math:`G`
   * - ``rad.kappa``
     - 1
     - Absorption coefficient :math:`\kappa`
   * - ``rad.emis``
     - 1
     - Emission :math:`\kappa B`

Spray (``PELE_USE_SPRAY``)
^^^^^^^^^^^^^^^^^^^^^^^^^^

Published when compiled with ``USE_PARTICLES=TRUE``. Spray derived quantities
are projected onto the AMR mesh grid via virtual particle interpolation.
Field names are assigned by the spray model and are mechanism-dependent;
the exact names can be discovered at runtime (see :ref:`sec:insitu::discovery`).
Typical fields include:

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``spray_num``
     - 1
     - Droplet number density (parcels per unit volume)
   * - ``spray_mass``
     - 1
     - Spray mass density
   * - ``spray_vol_frac``
     - 1
     - Liquid volume fraction
   * - ``spray_temp``
     - 1
     - Droplet temperature
   * - ``spray_density``
     - 1
     - Droplet material density
   * - ``d10``
     - 1
     - Arithmetic mean diameter
   * - ``d32``
     - 1
     - Sauter mean diameter
   * - ``spray_x_vel``, ``spray_y_vel``
     - 1 each
     - Droplet velocity components
   * - ``wall_film_hght``, ``wall_film_mass``
     - 1 each
     - Wall film height and mass (when wall film model is active)

Plasma (``PELE_USE_PLASMA``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Published when compiled with ``USE_PLASMA=TRUE``.

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``nE``
     - 1
     - Electron number density
   * - ``phiV``
     - 1
     - Electric potential
   * - ``DriftFlux_<ion>_X`` [``_Y``, ``_Z``]
     - NUM_IONS × SPACEDIM
     - Ion drift fluxes (when ``peleLM.do_extraEFdiags = 1``)

Embedded boundary (``AMREX_USE_EB``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Published when compiled with ``USE_EB=TRUE``.

.. list-table::
   :widths: 35 15 50
   :header-rows: 1

   * - Field name
     - Components
     - Description
   * - ``volFrac``
     - 1
     - EB volume fraction

.. _sec:insitu::discovery:

Discovering available field names at runtime
--------------------------------------------

Because some field names depend on the chemistry mechanism or physics model
(e.g. species names in ``rho.Y(<species>)`` and ``I_R(<species>)``, or soot
and spray model variable names), it is useful to have Ascent report exactly
which fields are present in the published mesh. This can be done by requesting
a nonexistent field in ``ascent_actions.yaml``: ::

    -
      action: "add_scenes"
      scenes:
        s1:
          plots:
            p1:
              type: "pseudocolor"
              field: "DISCOVER_FIELDS"
          renders:
            r1:
              image_prefix: "discover_%05d"
              image_width: 512
              image_height: 512

Ascent will print a line of the form: ::

    (s1/p1) unknown field 'DISCOVER_FIELDS' field names: 'RhoRT', 'Y(N2)',
    'avg_pressure', 'density', 'soot_fv', 'soot_N', 'rad.G', ...

listing every field available for that run. This is the canonical way to
obtain the exact field name strings for a given mechanism and physics
configuration before writing a production actions file.

.. _sec:insitu::examples:

Example actions files
---------------------

Pseudocolor temperature with AMR mesh overlay
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

    -
      action: "add_scenes"
      scenes:
        s1:
          image_prefix: "temp_mesh_%05d"
          plots:
            plt1:
              type: "pseudocolor"
              field: "temp"
            plt2:
              type: "mesh"

Multiple fields in a single run
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Multiple scenes are rendered simultaneously at each in-situ call with no
additional solver cost — the mesh is published once and all scenes consume
it. ::

    -
      action: "add_scenes"
      scenes:
        s1:
          image_prefix: "temp_%05d"
          plots:
            p1:
              type: "pseudocolor"
              field: "temp"
        s2:
          image_prefix: "pressure_%05d"
          plots:
            p1:
              type: "pseudocolor"
              field: "avg_pressure"
        s3:
          image_prefix: "vorticity_%05d"
          plots:
            p1:
              type: "pseudocolor"
              field: "mag_vort"

Reaction rate visualization (reacting cases)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Requires ``peleLM.plot_react = 1`` in the input file. ::

    -
      action: "add_scenes"
      scenes:
        s1:
          image_prefix: "heatrelease_%05d"
          plots:
            p1:
              type: "pseudocolor"
              field: "HeatRelease"
        s2:
          image_prefix: "IR_CH4_%05d"
          plots:
            p1:
              type: "pseudocolor"
              field: "I_R(CH4)"

.. _sec:insitu::example_run:

Example run
-----------

To run the ``HotBubble`` case with in-situ rendering every 50 steps,
suppressing plotfile and checkpoint output: ::

    mpirun -n 1 ./PeleLMeX2d.gnu.MPI.CUDA.ex input.2d-regt \
        amr.max_step=400 amr.plot_int=-1 amr.check_int=-1   \
        ascent.plot_int=50

A reference ``ascent_actions.yaml`` is provided in the ``HotBubble`` case
directory at ``Exec/RegTests/HotBubble/ascent_actions.yaml``.
