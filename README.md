# PeleLMeX

[![AMReX Badge](https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22)](https://amrex-codes.github.io/amrex/)
[![Exascale Computing Project](https://img.shields.io/badge/supported%20by-ECP-blue)](https://www.exascaleproject.org/research-project/combustion-pele/)
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-blue)](https://isocpp.org/)
[![Citing](https://joss.theoj.org/papers/10.21105/joss.05450/status.svg)](https://joss.theoj.org/papers/10.21105/joss.05450)
[![Archive](https://zenodo.org/badge/DOI/10.5281/zenodo.10056232.svg)](https://doi.org/10.5281/zenodo.10056232)

![CI](https://github.com/AMReX-Combustion/PeleLMeX/workflows/PeleLMeX-CI/badge.svg)

## Overview

*PeleLMeX* is a solver for high fidelity reactive flow simulations, namely direct numerical simulation (DNS) and large eddy simulation (LES).
The solver combines a low Mach number approach, adaptive mesh refinement (AMR), embedded boundary (EB) geometry treatment and high performance
computing (HPC) to provide
a flexible tool to address research questions on platforms ranging from small workstations to the world's largest GPU-accelerated supercomputers.
*PeleLMeX* has been used to study complex flame/turbulence interactions in [RCCI engines](https://gfm.aps.org/meetings/dfd-2022/63236765199e4c2c0873f9f6) and
[hydrogen combustion](https://www.sciencedirect.com/science/article/pii/S001021802300192X) or [the effect of sustainable aviation fuel
on gas turbine combustion](https://www.osti.gov/biblio/1995457).

*PeleLMeX* is part of the [Pele combustion Suite](https://amrex-combustion.github.io/).

## Documentation

![Documentation](https://github.com/AMReX-Combustion/PeleLMeX/workflows/PeleLMeX-Docs/badge.svg)

*PeleLMeX* is a non-subcycling version of [*PeleLM*](https://github.com/AMReX-Combustion/PeleLM) based on AMReX's
[AmrCore](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html) and borrowing from the incompressible
solver [incflo](https://github.com/AMReX-Codes/incflo). It solves of the multispecies reactive Navier-Stokes equations
in the low Mach number limit as described in the [documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/index.html).
It inherits most of *PeleLM* algorithmic features, but differs significantly in its implementation stemming from the non-subcycling approach.
*PeleLM* is no longer under active development; *PeleLMeX* should be used for simulations of low Mach number reacting flows and
[*PeleC*](https://github.com/AMReX-Combustion/PeleC) for simulations of flows with higher Mach numbers where compressibility effects are
significant.

A overview of *PeleLMeX* controls is provided in the [documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/LMeXControls.html).

### Core Algorithm

The *PeleLMeX* governing equations and core algorithms are described in:

https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#mathematical-background

https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#pelelmex-algorithm

### Tutorials

A set of self-contained tutorials describing more complex problems is also provided:

https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials.html

## Installation

### Requirements

The compilations of *PeleLMeX* requires a C++17 compatible compiler (GCC >= 8 or Clang >= 3.6) as
well as [CMake](https://cmake.org/download/) >= 3.23 for compiling the [SUNDIALS](https://github.com/LLNL/sundials) third party library.

Most of the examples provided hereafter and in the [tutorials](https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials.html)
will use MPI to run in parallel. Although not mandatory, it is advised to build *PeleLMeX* with MPI support from the get go if
more than a single core is available to you. Any of [mpich](https://www.mpich.org/downloads/) or
[open-mpi](https://www.open-mpi.org/software/ompi/v4.1/) is a suitable option if MPI is not already available on your platform.

Finally, when building with GPU support, CUDA >= 11 is required with NVIDIA GPUs and ROCm >= 5.2 is required with AMD GPUs.

### Download

The preferred method consists of cloning *PeleLMeX* and its submodules
([PelePhysics](https://github.com/AMReX-Combustion/PelePhysics),
[amrex](https://github.com/AMReX-Codes/amrex),
[AMReX-Hydro](https://github.com/AMReX-Fluids/AMReX-Hydro), and
[SUNDIALS](https://github.com/LLNL/sundials) using a recursive `git clone`:

```
git clone --recursive --shallow-submodules --single-branch https://github.com/AMReX-Combustion/PeleLMeX.git
```

The `--shallow-submodules` and `--single-branch` flags are recommended for most users as they substantially reduce the size of the download by
skipping extraneous parts of the git history. Developers may wish to omit these flags in order download the complete git history of PeleLMeX
and its submodules, though standard `git` commands may also be used after a shallow clone to obtain the skipped portions if needed.

Alternatively, you can use a separate `git clone` of each of the submodules.
The default location for *PeleLMeX* dependencies is the `Submodules` folder but you optionally
setup the following environment variables (e.g. using bash) to any other location:

```
export PELE_HOME=<path_to_PeleLMeX>
export AMREX_HYDRO_HOME=${PELE_HOME}/Submodules/AMReX-Hydro
export PELE_PHYSICS_HOME=${PELE_HOME}/Submodules/PelePhysics
export AMREX_HOME=${PELE_PHYSICS_HOME}/Submodules/amrex
export SUNDIALS_HOME=${PELE_PHYSICS_HOME}/Submodules/sundials
```

### Compilation

Both GNUmake and CMake can be used to build *PeleLMeX* executables. GNUmake is the preferred choice for single executables when running production simulations. While CMake is the preferred method for automatically building and testing most available executables.
The code handling the initial condition and boundary conditions is unique to each case,
and subfolders in the `Exec` directory provide a number of examples.

For instance, to compile the executable for the case of a rising hot bubble,
move into the `HotBubble` folder:

```
cd PeleLMeX/Exec/RegTests/HotBubble
```

If this is a clean install, you will need to make the third party libraries with: `make TPL` (note: if on macOS, you might need to specify `COMP=llvm` in the `make` statements).

Finally, make with: `make -j`, or if on macOS: `make -j COMP=llvm`. To clean the installation, use either `make clean` or `make realclean`. If running into compile errors after changing compile time options in PeleLMeX (e.g., the chemical mechanism), the first thing to try is to clean your build by running `make TPLrealclean && make realclean`, then try to rebuild the third party libraries and PeleLMeX with `make TPL && make -j`. See the [Tutorial](https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials_HotBubble.html) for this case for instructions on how to compile with different options (for example, to compile without MPI support or to compile for GPUs) and how to run the code once compiled.

To compile and test using CMake, refer to the example `cmake.sh` script in the `Build` directory, or reference the GitHub Actions workflows in the `.github/workflows` directory.

## Getting help, contributing

Do you have a question ? Found an issue ? Please use the [GitHub Discussions](https://github.com/AMReX-Combustion/PeleLMeX/discussions) to engage
with the development team or open a new [GitHub issue](https://github.com/AMReX-Combustion/PeleLMeX/issues) to report a bug. The development team
also encourages users to take an active role in respectfully answering each other's questions in these spaces. When reporting a bug, it is helpful
to provide as much detail as possible, including a case description and the major compile and runtime options being used. Though not required,
it is most effective to create a fork of this repository and share a branch of that fork with a case that minimally reproduces the error.

New contributions to *PeleLMeX* are welcome ! Contributing Guidelines are provided in [CONTRIBUTING.md](CONTRIBUTING.md).

## Acknowledgment

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.

## Citation

To cite PeleLMeX, please use [![Citing](https://joss.theoj.org/papers/10.21105/joss.05450/status.svg)](https://joss.theoj.org/papers/10.21105/joss.05450)

```
@article{PeleLMeX_JOSS,
  doi = {10.21105/joss.05450},
  url = {https://doi.org/10.21105/joss.05450},
  year = {2023},
  month = october,
  publisher = {The Open Journal},
  volume = {8},
  number = {90},
  pages = {5450},
  author = {Lucas Esclapez and Marc Day and John Bell and Anne Felden and Candace Gilet and Ray Grout and Marc {Henry de Frahan} and Emmanuel Motheau and Andrew Nonaka and Landon Owen and Bruce Perry and Jon Rood and Nicolas Wimer and Weiqun Zhang},
  journal = {Journal of Open Source Software}
}
```
