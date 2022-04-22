TOP = ../../../..
PELELMEX_HOME      ?= ${TOP}/PeleLMeX
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
# AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex_original
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro
# HYPRE_HOME         ?= ${PELELMEX_HOME}/Submodules/hypre
HYPRE_HOME          = /lustre/scratch/bsouzas/PeleLMeX/Submodules/hypre/install

# AMReX
DIM             = 3
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = TRUE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE
USE_EB          = TRUE
USE_HYPRE       = TRUE
CEXE_headers 	+= EBUserDefined.H

# PelePhysics
Chemistry_Model = air
Eos_Model       = Fuego
Transport_Model = Simple

ifeq ($(USE_CUDA),TRUE)
DEFINES += -DAMREX_GPUS_PER_NODE=2
endif

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
