# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_EB = FALSE
# input.2d_ShearLayer_beta2 needs this TRUE: at cosh^2(2.4) = 30.9 the
# default AMReX MLMG projections fail immediately.  Every other input in
# this directory runs fine with it FALSE.
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

# PelePhysics
Chemistry_Model = air
Eos_Model = Fuego
Transport_Model = Constant

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX
