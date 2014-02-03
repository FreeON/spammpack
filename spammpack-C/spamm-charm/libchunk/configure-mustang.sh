#!/bin/bash

module purge
module load intel
module load mkl

MKL_CPPFLAGS="-I${MKLROOT}/include"
#MKL_CFLAGS="-mkl=sequential"
MKL_CFLAGS="-mkl=parallel"

./configure \
  --enable-block-multiply=blas \
  --disable-assert \
  CC=icc \
  CPPFLAGS=${MKL_CPPFLAGS} \
  CFLAGS=${MKL_CFLAGS}
