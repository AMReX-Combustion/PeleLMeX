TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro


# AMReX
DIM             = 3
DEBUG           = TRUE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE
USE_SYCL        = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
Chemistry_Model = dodecane_lu
Eos_Model		= Fuego
Transport_Model = Simple

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
