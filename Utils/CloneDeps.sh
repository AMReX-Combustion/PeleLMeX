#!/usr/bin/env bash

echo "Getting PeleLMeX dependencies - tests ... "
export PELELM_HOME=${PWD}/..
mkdir build
git clone https://github.com/AMReX-Codes/amrex.git build/amrex
export AMREX_HOME=${PWD}/build/amrex
git clone -b development https://github.com/AMReX-Combustion/PelePhysics.git build/PelePhysics
export PELE_PHYSICS_HOME=${PWD}/build/PelePhysics
git clone https://github.com/AMReX-Codes/AMReX-Hydro.git build/AMReX-Hydro
export AMREX_HYDRO_HOME=${PWD}/build/AMReX-Hydro
git clone https://github.com/AMReX-Combustion/PeleMP.git build/PeleMP
export PELEMP_HOME=${PWD}/build/PeleMP
git clone https://github.com/LLNL/sundials.git build/sundials
export SUNDIALS_HOME=${PWD}/build/sundials
