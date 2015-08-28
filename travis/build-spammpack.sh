#!/bin/bash -x

cmake . \
  -DCMAKE_INSTALL_PREFIX=${PWD}/install \
  -DCMAKE_C_COMPILER=${CC:=gcc-${GCC_VERSION}} \
  -DCMAKE_Fortran_COMPILER=${FC:=gfortran-${GCC_VERSION}} \
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
  -DCMAKE_VERBOSE_MAKEFILE=yes \
  -DSPAMM_TESTING=yes \
  || exit
make || exit
make install || exit
ctest --verbose || exit
