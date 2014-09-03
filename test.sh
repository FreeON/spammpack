#!/bin/bash

FC=gfortran CC=gcc CMAKE_BUILD_TYPE=Release ./tester.sh || exit
FC=gfortran CC=gcc CMAKE_BUILD_TYPE=Debug ./tester.sh || exit

. /opt/intel/bin/compilervars.sh intel64
FC=ifort CC=icc CMAKE_BUILD_TYPE=Release ./tester.sh || exit
FC=ifort CC=icc CMAKE_BUILD_TYPE=Debug ./tester.sh || exit
