TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
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

# PeleLMeX
USE_EFIELD		 = TRUE

# PelePhysics
Chemistry_Model = IonizedAir
Reactions_dir	 = Fuego
Eos_dir			 = Fuego
Transport_dir   = Simple

USE_SUNDIALS_PP = TRUE

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
