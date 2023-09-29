# PeleLMeX

[![AMReX Badge](https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22)](https://amrex-codes.github.io/amrex/)
[![Exascale Computing Project](https://img.shields.io/badge/supported%20by-ECP-blue)](https://www.exascaleproject.org/research-project/combustion-pele/)
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-blue)](https://isocpp.org/)
[![JOSS](https://joss.theoj.org/papers/6142eb838783b07fce14450fefe21e07/status.svg)](https://joss.theoj.org/papers/6142eb838783b07fce14450fefe21e07)

![CUDA build](https://github.com/AMReX-Combustion/PeleLMeX/workflows/PeleLMeX_Cuda/badge.svg)
![HIP build](https://github.com/AMReX-Combustion/PeleLMeX/workflows/PeleLMeX_Hip/badge.svg)
![SYCL build](https://github.com/AMReX-Combustion/PeleLMeX/workflows/PeleLMeX_Intel/badge.svg)

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
([amrex](https://github.com/AMReX-Codes/amrex), [AMReX-Hydro](https://github.com/AMReX-Codes/AMReX-Hydro), [PelePhysics](https://github.com/AMReX-Combustion/PelePhysics), [PeleMP](https://github.com/AMReX-Combustion/PeleMP)), and [SUNDIALS](https://github.com/LLNL/sundials) using a recursive `git clone`:

```
git clone --recursive https://github.com/AMReX-Combustion/PeleLMeX.git
```

Alternatively, you can use a separate `git clone` of each of the submodules.
The default location for *PeleLMeX* dependencies is the `Submodules` folder but you optionally
setup the following environment variables (e.g. using bash) to any other location:

```
export PELELMEX_HOME=<path_to_PeleLMeX>
export AMREX_HOME=${PELELMEX_HOME}/Submodules/amrex
export AMREX_HYDRO_HOME=${PELELMEX_HOME}/Submodules/AMReX-Hydro
export PELE_PHYSICS_HOME=${PELELMEX_HOME}/Submodules/PelePhysics
export PELEMP_HOME=${PELELMEX_HOME}/Submodules/PeleMP
export SUNDIALS_HOME=${PELELMEX_HOME}/Submodules/sundials
```

### Compilation

Both GNUmake and CMake can be used to build a *PeleLMeX* executable, but GNUmake is the preferred choice.
The code handling the initial condition and boundary conditions is unique to each case,
and subfolders in the `Exec` directory provide a number of examples.

For instance, to compile the executable for the case of a rising hot bubble,
move into the `HotBubble` folder:

```
cd PeleLMeX/Exec/RegTests/HotBubble
```

If this is a clean install, you will need to make the third party libraries with: `make TPL` (note: if on macOS, you might need to specify `COMP=llvm` in the `make` statements).

Finally, make with: `make -j`, or if on macOS: `make -j COMP=llvm`. To clean the installation, use either `make clean` or `make realclean`. If running into compile errors after changing compile time options in PeleLMeX (e.g., the chemical mechanism), the first thing to try is to clean your build by running `make TPLrealclean && make realclean`, then try to rebuild the third party libraries and PeleLMeX with `make TPL && make -j`. See the [Tutorial](https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials_HotBubble.html) for this case for instructions on how to compile with different options (for example, to compile without MPI support or to compile for GPUs) and how to run the code once compiled.

When using CMake (cmake version >= 3.23 is required), first configure CMake directly in *PeleLMeX* root folder:

```
cmake -S . -B buildHotBubble -DPELELMEX_MPI=ON -DPELELMEX_CASE=HotBubble
```

Then build the executable:

```
cmake --build buildHotBubble --parallel 4
```

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
