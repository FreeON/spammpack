#!/bin/bash

git clean -df

cmake . \
  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:-Debug} \
  -DDEBUG_LEVEL=2 \
  -DCMAKE_VERBOSE_MAKEFILE=yes \
  -DCMAKE_C_COMPILER=${CC:-gcc} \
  -DCMAKE_Fortran_COMPILER=${FC:-gfortran} || exit