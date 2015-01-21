#!/bin/bash

mkdir -p build
cd build

cmake .. \
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Debug} \
  -DSPAMM_DEBUG_LEVEL=${SPAMM_DEBUG_LEVEL:-2} \
  -DCMAKE_VERBOSE_MAKEFILE=yes \
  -DCMAKE_C_COMPILER=${CC:-gcc} \
  -DCMAKE_Fortran_COMPILER=${FC:-gfortran} || exit

echo "The sources are configured in the build directory."
