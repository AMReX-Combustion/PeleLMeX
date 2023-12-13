# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
TINY_PROFILE    = FALSE

# Compilation
COMP            = llvm
USE_MPI         = FALSE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE
USE_SYCL        = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
Chemistry_Model = SootReaction
Eos_Model	= Fuego
Transport_Model = Simple

USE_PELERAD = TRUE
USE_SOOT = TRUE
NUM_SOOT_MOMENTS = 3

PELELMEX_HOME ?= ../../..
include $(PELELMEX_HOME)/Exec/Make.PeleLMeX
