# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = TRUE

# Compilation
COMP            = llvm
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD      = FALSE

# PelePhysics
Chemistry_Model = SootReaction
Eos_Model	= Fuego
Transport_Model = Simple

USE_SOOT = TRUE
NUM_SOOT_MOMENTS = 3

USE_PELERAD = TRUE
ifeq ($(USE_PELERAD), TRUE)
ifeq ($(USE_HIP), TRUE)
LIBRARIES += -lstdc++fs
endif
endif

include $(PELELMEX_HOME)/Exec/Make.PeleLMeX
