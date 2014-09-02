#!/bin/bash

git clean -df

cmake . -DCMAKE_C_COMPILER=gcc -DCMAKE_Fortran_COMPILER=gfortran || exit
make || exit
make test || exit

. /opt/intel/bin/compilervars.sh intel64

git clean -df

cmake . -DCMAKE_C_COMPILER=icc -DCMAKE_Fortran_COMPILER=ifort || exit
make || exit
make test || exit
