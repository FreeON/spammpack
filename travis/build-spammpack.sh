#!/bin/bash -x

cmake . \
  -DCMAKE_INSTALL_PREFIX=${PWD}/install \
  -DCMAKE_Fortran_COMPILER=gfortran-${GCC_VERSION} \
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
  -DCMAKE_VERBOSE_MAKEFILE=yes || exit
make || exit
make install || exit
ctest --verbose || exit
