TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELM_HOME        ?= ${TOP}/PeleLM
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics


# AMReX
DIM             = 2
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

# PeleLM

# PelePhysics
Chemistry_Model = air
Reactions_dir	 = Null
Eos_dir			 = Fuego
Transport_dir   = Simple

include $(PELELM_HOME)/Utils/Make.PeleLM
