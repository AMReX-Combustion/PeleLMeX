cmake_minimum_required(VERSION 3.14)

project( PeleLMeX
         DESCRIPTION "An AMR Low Mach Number Reactive Flow Simulation Code Without Level Sub-cycling"
         HOMEPAGE_URL "https://github.com/AMReX-Combustion/PeleLMeX"
         LANGUAGES CXX
       )

message(STATUS "CMake version: ${CMAKE_VERSION}")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Utils/CMake")

#
# Base options
#
set(PELELMEX_DIM             "3"      CACHE STRING "Number of physical dimensions")
set(PELELMEX_PRECISION       "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")
set(PELELMEX_CHEMISTRY_MODEL drm19    CACHE STRING "Chemical mechanism")
set(PELELMEX_EOS_MODEL       Fuego    CACHE STRING "Equation of state model")
set(PELELMEX_TRANSPORT_MODEL Simple   CACHE STRING "Transport model")

#
# Physics options
#
option(PELELMEX_USE_EB        "Enable Embedded Boundary"     OFF)
option(PELELMEX_USE_EFIELD    "Enable Electric field module" OFF)
option(PELELMEX_USE_PARTICLES "Enable Spray module"          OFF)
option(PELELMEX_USE_SOOR      "Enable Soot module"           OFF)

#
# Misc options
#
option(PELELMEX_ENABLE_HYPRE        "Enable Offload to Hypre"    OFF)
option(PELELMEX_ENABLE_TINY_PROFILE "Enable AMReX tiny profiler" OFF)
option(PELELMEX_ENABLE_HDF5         "Enable HDF5 PlotFiles"      OFF)

#
# HPC options
#
option(PELELMEX_MPI    "Enable MPI"    OFF)
option(PELELMEX_OPENMP "Enable OpenMP" OFF)
option(PELELMEX_CUDA   "Enable CUDA"   OFF)
option(PELELMEX_HIP    "Enable HIP"    OFF)
option(PELELMEX_SYCL   "Enable SyCL"   OFF)

# C++ Options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Checks
if(PELELMEX_DIM EQUAL 1)
    message(FATAL_ERROR "PeleLMeX does not support 1D.")
endif()

########################### PeleLMeX ##################################
add_executable(pelelmex "")
set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)
set(EXEC_DIR ${CMAKE_SOURCE_DIR}/Exec/RegTests/FlameSheet)
add_subdirectory(${SRC_DIR})
add_subdirectory(${EXEC_DIR})
set(PELEPHYSICS_LIB PelePhysics-${PELELMEX_CHEMISTRY_MODEL}-${PELELMEX_EOS_MODEL}-${PELELMEX_TRANSPORT_MODEL})
include(SetRpath)
target_link_libraries(pelelmex
    PUBLIC
        AMReX-Hydro::amrex_hydro_api
        AMReX::amrex
        SUNDIALS::cvode
        SUNDIALS::arkode
        PelePhysics::${PELEPHYSICS_LIB}_api)

########################### AMReX #####################################

message(STATUS "AMReX Configure Section")
set(AMREX_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/amrex")
include(SetAmrexOptions)
list(APPEND CMAKE_MODULE_PATH "${AMREX_SUBMOD_LOCATION}/Tools/CMake")
add_subdirectory(${AMREX_SUBMOD_LOCATION})
include(AMReXBuildInfo)
generate_buildinfo(pelelmex ${CMAKE_SOURCE_DIR})
target_include_directories(pelelmex SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

########################### AMReX-Hydro ###############################

message(STATUS "AMReX-Hydro Configure Section")
set(AMREX_HYDRO_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/AMReX-Hydro")
include(SetAmrexHydroOptions)
add_subdirectory(${AMREX_HYDRO_SUBMOD_LOCATION})

########################### SUNDIALS ##################################

message(STATUS "Sundials Configure Section")
set(SUNDIALS_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/sundials")
include(SetSundialsOptions)
add_subdirectory(${SUNDIALS_SUBMOD_LOCATION})

########################### PelePhysics ###############################

message(STATUS "PelePhysics Configure Section")
set(PELE_PHYSICS_SUBMOD_LOCATION "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics")
include(SetPelePhysicsOptions)
add_subdirectory(${PELE_PHYSICS_SUBMOD_LOCATION})
