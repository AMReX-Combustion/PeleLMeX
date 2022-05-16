TOP = ../../../..
AMREX_HOME         ?= ${TOP}/amrex
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
PELE_PHYSICS_HOME  ?= ${TOP}/PelePhysics
AMREX_HYDRO_HOME   ?= ${TOP}/AMReX-Hydro
HYPRE_HOME         = /ccs/home/lesclapez/hypre/install


# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE
USE_EB			 = TRUE
USE_HYPRE		 = TRUE

# Compilation
COMP            = clang
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = TRUE

# PeleLMeX
CEXE_headers 	+= EBUserDefined.H

# PelePhysics
Chemistry_Model = dodmethair_4sp
Eos_Model		 = Fuego
Transport_Model = Simple

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
