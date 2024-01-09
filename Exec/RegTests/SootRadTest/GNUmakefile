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

USE_RADIATION = TRUE
USE_SOOT = TRUE
NUM_SOOT_MOMENTS = 3

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX
