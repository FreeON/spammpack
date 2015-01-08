#!/bin/bash

module purge
module load intel
module load mkl

MKL_CPPFLAGS="-I${MKLROOT}/include"

if [[ $1 == "parallel" ]]; then
  MKL_CFLAGS="-mkl=parallel"
else
  MKL_CFLAGS="-mkl=sequential"
fi

./configure \
  --enable-block-multiply=blas \
  --disable-assert \
  CC=icc \
  CPPFLAGS=${MKL_CPPFLAGS} \
  CFLAGS=${MKL_CFLAGS}
