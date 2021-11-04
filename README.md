# PeleLMeX

[![AMReX Badge](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io/amrex/)

## Overview

*PeleLMeX* is a non-subcycling version of [*PeleLM*](https://github.com/AMReX-Combustion/PeleLM) based on AMReX's [AmrCore](https://amrex-codes.github.io/amrex/docs_html/AmrCore.html) and borrowing from the incompressible solver [incflo](https://github.com/AMReX-Codes/incflo).

## Installation

First, clone the prerequisite AMReX-Codes and AMReX-Combustion repositories:

amrex: `git clone https://github.com/AMReX-Codes/amrex.git`

AMReX-Hydro: `git clone https://github.com/AMReX-Codes/AMReX-Hydro.git`

PelePhysics: `git clone https://github.com/AMReX-Combustion/PelePhysics.git`

Next, clone the PeleLMeX repository:

`git clone https://github.com/AMReX-Combustion/PeleLMeX.git`

`cd PeleLMeX`

Ensure the following paths are set somewhere (e.g., in `.bashrc`) and point to the recently cloned repos:
```
AMREX_HOME=${Path_to_amrex}
AMREX_HYDRO_HOME=${Path_to_AMReX-Hydro}
PELE_PHYSICS_HOME=${Path_to_PelePhysics}
PELELMEX_HOME=${Path_to_PeleLMeX}
```

Move to a sample case to compile (such as `HotBubble`):

`cd Exec/RegTests/HotBubble`

If this is a clean install, you will need to make the third party libraries with: `make TPL` (note: if on macOS, you might need to specify `COMP=llvm` in the `make` statements).

Finally, make with: `make -j`, or if on macOS: `make -j COMP=llvm`.

To clean the installation, use either `make clean` or `make realclean`.