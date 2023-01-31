TOP = ../../..
PELELMEX_HOME      ?= ${TOP}
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro

# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE
USE_EB			= TRUE
USE_HYPRE		= FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX

# PelePhysics
Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Constant

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
