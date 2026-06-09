# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_EB = FALSE
# HYPRE / BoomerAMG is the recommended MAC bottom solver at beta >= 2.
# Default to TRUE because the shipped input files (beta = 2 2) require
# it.  Set to FALSE if hypre is not available; in that case also drop
# TanhStretchMap.beta to 1.0 1.0 in the input file and remove the
# mac_proj.use_mlhypre / nodal_proj.bottom_solver lines.
USE_HYPRE = TRUE

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

# PelePhysics (incompressible; chemistry/EOS unused but needed to link)
Chemistry_Model = air
Eos_Model = Fuego
Transport_Model = Simple

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX
