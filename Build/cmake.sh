#!/bin/bash

# Example CMake config script for running test suite on a MacOS laptop with OpenMPI:
#export PATH=/opt/homebrew/Cellar/llvm/16.0.6/bin:${PATH}
cmake -DCMAKE_INSTALL_PREFIX:PATH=./install \
      -DCMAKE_CXX_COMPILER:STRING=clang++ \
      -DCMAKE_C_COMPILER:STRING=clang \
      -DMPIEXEC_PREFLAGS:STRING=--oversubscribe \
      -DCMAKE_BUILD_TYPE:STRING=Release \
      -DPELELMEX_DIM:STRING=3 \
      -DPELELMEX_ENABLE_MPI:BOOL=OFF \
      -DPELELMEX_ENABLE_FCOMPARE:BOOL=OFF \
      -DPELELMEX_ENABLE_FCOMPARE_FOR_TESTS:BOOL=OFF \
      -DPELELMEX_ENABLE_MASA:BOOL=OFF \
      -DPELELMEX_ENABLE_CPPCHECK:BOOL=OFF \
      -DPELELMEX_ENABLE_CLANG_TIDY:BOOL=OFF \
      -DPELELMEX_ENABLE_CUDA:BOOL=OFF \
      -DAMReX_CUDA_ARCH=Volta \
      -DPYTHON_EXECUTABLE=$(which python3) \
      -DPELELMEX_PRECISION:STRING=DOUBLE \
      ..
#make
cmake --build . --parallel $(sysctl -n hw.ncpu) #&> output.txt
#ctest -j $(sysctl -n hw.ncpu)
