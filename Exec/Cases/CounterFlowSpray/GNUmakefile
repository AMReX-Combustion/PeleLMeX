TOP = ../../..
PELELMEX_HOME      ?= ${TOP}
AMREX_HOME         ?= ${PELELMEX_HOME}/Submodules/amrex
PELE_PHYSICS_HOME  ?= ${PELELMEX_HOME}/Submodules/PelePhysics
AMREX_HYDRO_HOME   ?= ${PELELMEX_HOME}/Submodules/AMReX-Hydro
PELEMP_HOME        ?= ${PELELMEX_HOME}/Submodules/PeleMP

# AMReX
DIM             = 2
DEBUG           = FALSE
PRECISION       = DOUBLE
VERBOSE         = FALSE
TINY_PROFILE    = FALSE

# Compilation
COMP            = gnu
USE_MPI         = TRUE
USE_OMP         = FALSE
USE_CUDA        = FALSE
USE_HIP         = FALSE

# PeleLMeX
USE_EFIELD		 = FALSE

# PeleMP
USE_PARTICLES    = TRUE
SPRAY_FUEL_NUM = 1
CEXE_sources += SprayParticlesInitInsert.cpp

# PelePhysics
Chemistry_Model = heptane_fc
Eos_Model	    = Fuego
Transport_Model = Simple

include $(PELELMEX_HOME)/Utils/Make.PeleLMeX
