TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro


# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE
USE_EB 			= TRUE
USE_HYPRE		= FALSE

# Compilation
COMP            = llvm
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE
USE_SYCL        = FALSE

# PeleLMeX

# PelePhysics
Chemistry_Model = drm19
Eos_Model		= Fuego
Transport_Model = Simple

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
