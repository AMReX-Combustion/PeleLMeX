# AMReX
DIM = 2
COMP = gnu
PRECISION = DOUBLE
USE_EB = FALSE
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

# PelePhysics
USE_PLASMA = FALSE
USE_PARTICLES = TRUE
SPRAY_GCM = TRUE
SPRAY_FUEL_NUM = 2
ifeq ($(SPRAY_GCM), TRUE)
  SpraySuffix = .$(SPRAY_FUEL_NUM)SprayFuel.SprayGCM
else
  SpraySuffix = .$(SPRAY_FUEL_NUM)SprayFuel.SprayMP
endif

USE_MANIFOLD = FALSE
ifeq ($(USE_MANIFOLD), TRUE)
  Eos_Model       = Manifold
  Chemistry_Model = Null
  Transport_Model = Manifold
  Manifold_Dim = $(SPRAY_FUEL_NUM)
  Manifold_Type = Table
  USERSuffix = .Manifold$(Manifold_Dim)d$(Manifold_Type)$(SpraySuffix)
else
  SPRAY_FUEL_NUM = 2
  Chemistry_Model = liquid_fuels_nonreacting
  Eos_Model = Fuego
  Transport_Model = Simple
  USERSuffix = .$(Chemistry_Model)$(SpraySuffix)
endif

PELE_HOME ?= ../../..
include $(PELE_HOME)/Exec/Make.PeleLMeX
