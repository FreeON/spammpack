#!/bin/bash

. /opt/intel/bin/compilervars.sh intel64

FC=ifort CC=icc CMAKE_BUILD_TYPE=Debug ./configure.sh || exit
make || exit
make test || exit

FC=ifort CC=icc CMAKE_BUILD_TYPE=Release ./configure.sh || exit
make || exit
make test || exit

FC=gfortran CC=gcc CMAKE_BUILD_TYPE=Debug ./configure.sh || exit
make || exit
make test || exit

FC=gfortran CC=gcc CMAKE_BUILD_TYPE=Release ./configure.sh || exit
make || exit
make test || exit
