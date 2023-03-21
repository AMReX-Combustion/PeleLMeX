# PeleLMeX

[![AMReX Badge](https://img.shields.io/static/v1?label=%22powered%20by%22&message=%22AMReX%22&color=%22blue%22)](https://amrex-codes.github.io/amrex/)
[![Exascale Computing Project](https://img.shields.io/badge/supported%20by-ECP-blue)](https://www.exascaleproject.org/research-project/combustion-pele/)
[![Language: C++17](https://img.shields.io/badge/language-C%2B%2B17-blue)](https://isocpp.org/)

## Overview

![CUDA build](https://github.com/AMReX-Combustion/PeleLMeX/workflows/cuda/badge.svg)
![HIP build](https://github.com/AMReX-Combustion/PeleLMeX/workflows/hip/badge.svg)
![SYCL build](https://github.com/AMReX-Combustion/PeleLMeX/workflows/intel/badge.svg)

*PeleLMeX* is a non-subcycling version of [*PeleLM*](https://github.com/AMReX-Combustion/PeleLM) based on AMReX's [AmrCore](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html) and borrowing from the incompressible solver [incflo](https://github.com/AMReX-Codes/incflo). It is designed to run efficiently on small workstations as well as the largest ExaScale platforms currently available.

*PeleLMeX* is part of the [Pele combustion Suite](https://amrex-combustion.github.io/).

## Documentation

![Documentation](https://github.com/AMReX-Combustion/PeleLMeX/workflows/docs/badge.svg)

*PeleLMeX* solves of the multispecies reactive Navier-Stokes equations in the low Mach number limit as described in the [documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/index.html). It inherits most of *PeleLM* algorithmic features, but differs significantly in its implementation stemming from the non-subcycling approach.

A overview of *PeleLMeX* controls is provided in the [documentation](https://amrex-combustion.github.io/PeleLMeX/manual/html/LMeXControls.html).

### Core Algorithm

The *PeleLMeX* governing equations and core algorithms are described in:

https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#mathematical-background

https://amrex-combustion.github.io/PeleLMeX/manual/html/Model.html#pelelmex-algorithm

### Tutorials

A set of self-contained tutorials describing more complex problems is also provided:

https://amrex-combustion.github.io/PeleLMeX/manual/html/Tutorials.html

## Installation

### Download

The prefered method consist in cloning the *PeleLMeX* and its submodules using a recursive `git clone`:
([amrex](https://github.com/AMReX-Codes/amrex), [AMReX-Hydro](https://github.com/AMReX-Codes/AMReX-Hydro) and [PelePhysics](https://github.com/AMReX-Combustion/PelePhysics), [PeleMP](https://github.com/AMReX-Combustion/PeleMP)) using:

```
git clone --depth 1 --recursive https://github.com/AMReX-Combustion/PeleLMeX.git
```

Alternatively, you can use separate `git clone` of the each of the submodules.
The default location for *PeleLMeX* dependencies is the `Submodule` folder but you optionnally
setup the following environment variables (e.g. using bash) to an other location:

```
export PELELMEX_HOME=<path_to_PeleLMeX>
export AMREX_HOME=${PELELMEX_HOME}/Submodules/amrex
export AMREX_HYDRO_HOME=${PELELMEX_HOME}/Submodules/AMReX-Hydro
export PELE_PHYSICS_HOME=${PELELMEX_HOME}/Submodules/PelePhysics
export PELEMP_HOME=${PELELMEX_HOME}/Submodules/PeleMP
```

### Compilation

Both GNUmake and CMake can be used to build a *PeleLMeX* executable, but GNUmake is the prefered choice.
The code handling the initial condition and boundary conditions is unique to each case,
and subfolders in the `Exec` directory provide a number of example.

For instance, to compile the executable for the case of a rising hot bubble,
move into the `HotBubble` folder:

```
cd PeleLMeX/Exec/RegTest/HotBubble
```

If this is a clean install, you will need to make the third party libraries with: `make TPL` (note: if on macOS, you might need to specify `COMP=llvm` in the `make` statements).

Finally, make with: `make -j`, or if on macOS: `make -j COMP=llvm`. To clean the installation, use either `make clean` or `make realclean`.

When using CMake (cmake version >= 3.23 is required), first configure CMake directly in *PeleLMeX* root folder:

```
cmake -S . -B buildHotBubble -DPELELMEX_MPI=ON -DPELELMEX_CASE=HotBubble
```

Then build the executable:

```
cmake --build buildHoBubble --parallel 4
```

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
