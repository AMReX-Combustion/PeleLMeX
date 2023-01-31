# PeleLMeX

[![AMReX Badge](https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22)](https://amrex-codes.github.io/amrex/)
[![Exascale Computing Project](https://img.shields.io/badge/supported%20by-ECP-blue)](https://www.exascaleproject.org/research-project/combustion-pele/)
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-blue)](https://isocpp.org/)

## Overview

*PeleLMeX* is a non-subcycling version of [*PeleLM*](https://github.com/AMReX-Combustion/PeleLM) based on AMReX's [AmrCore](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html) and borrowing from the incompressible solver [incflo](https://github.com/AMReX-Codes/incflo).

*PeleLMeX* is part of the [Pele combustion Suite](https://amrex-combustion.github.io/).

## Documentation

*PeleLMeX* solves of the multispecies reactive Navier-Stokes equations in the low Mach number limit as described in the [documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/index.html). It inherits most of *PeleLM* algorithmic features, but differs significantly in its implementation stemming from the non-subcycling approach.

A overview of PeleLMeX controls is provided in the [documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/LMeXControls.html).

### Core Algorithm

The *PeleLMeX* governing equations and core algorithms are described in:

https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#mathematical-background

https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#pelelmex-algorithm

### Tutorials

A set of self-contained tutorials describing more complex problems is also provided:

https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials.html

## Installation

### Using git submodules

You can clone the PeleLMeX and tested versions of its submodules
([amrex](https://github.com/AMReX-Codes/amrex), [AMReX-Hydro](https://github.com/AMReX-Codes/AMReX-Hydro) and [PelePhysics](https://github.com/AMReX-Combustion/PelePhysics)) using:

```
git clone --recursive https://github.com/AMReX-Combustion/PeleLMeX.git
```

You can optionally setup the following environment variables (e.g. using bash):

```
export PELELMEX_HOME=<path_to_PeleLMeX>
export AMREX_HOME=${PELELMEX_HOME}/Submodules/amrex
export AMREX_HYDRO_HOME=${PELELMEX_HOME}/Submodules/AMReX-Hydro
export PELE_PHYSICS_HOME=${PELELMEX_HOME}/Submodules/PelePhysics
```

If you do not set these paths as environment variables, they will be assumed as specified in the `GNUmakefile` for the case you are compiling. If compiling in the default location for each case, no modifications are necessary.

Then, move into one of the available examples, such as `HotBubble`:

```
cd PeleLMeX/Exec/RegTest/HotBubble
```

If this is a clean install, you will need to make the third party libraries with: `make TPL` (note: if on macOS, you might need to specify `COMP=llvm` in the `make` statements).

Finally, make with: `make -j`, or if on macOS: `make -j COMP=llvm`.

### Using separate git clone

First, clone the source code and its prerequisite ([amrex](https://github.com/AMReX-Codes/amrex), [AMReX-Hydro](https://github.com/AMReX-Codes/AMReX-Hydro) and [PelePhysics](https://github.com/AMReX-Combustion/PelePhysics)) using:

```
git clone https://github.com/AMReX-Codes/amrex.git
git clone https://github.com/AMReX-Codes/AMReX-Hydro.git
git clone https://github.com/AMReX-Combustion/PelePhysics.git
git clone https://github.com/AMReX-Combustion/PeleLMeX.git
```

Ensure the following paths are set somewhere (e.g., in `.bashrc`) and point to the recently cloned repos:
```
export AMREX_HOME=<path_to_amrex>
export AMREX_HYDRO_HOME=<path_to_AMReX-Hydro>
export PELE_PHYSICS_HOME=<path_to_PelePhysics>
export PELELMEX_HOME=<path_to_PeleLMeX>
```

Move to a sample case to compile, such as `HotBubble`:

```
cd PeleLMeX/Exec/RegTests/HotBubble
```

If this is a clean install, you will need to make the third party libraries with: `make TPL` (note: if on macOS, you might need to specify `COMP=llvm` in the `make` statements).

Finally, make with: `make -j`, or if on macOS: `make -j COMP=llvm`.

To clean the installation, use either `make clean` or `make realclean`.

## Contributing

New contributions to *PeleLMeX* are welcome ! Contributing Guidelines are provided in [CONTRIBUTING.md](CONTRIBUTING.md).

## Acknowledgment

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
