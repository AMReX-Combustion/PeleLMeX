# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_EB = TRUE
USE_HYPRE = FALSE

# Profiling
PROFILE = FALSE
TINY_PROFILE = FALSE
COMM_PROFILE = FALSE
TRACE_PROFILE = FALSE
MEM_PROFILE = FALSE
USE_GPROF = FALSE

# Performance
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = FALSE
USE_HIP = FALSE
USE_SYCL = FALSE

# Debugging
DEBUG = FALSE
FSANITIZER = FALSE
THREAD_SANITIZER = FALSE

# PeleLMeX
CEXE_headers += EBUserDefined.H
CEXE_sources += PeleLMeX_PatchFlowVariables.cpp

# PelePhysics
Chemistry_Model = drm19
Eos_Model = Fuego
Transport_Model = Simple

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX
