#!/bin/bash

module purge
module load intel
module load mkl

./configure \
  --enable-block-multiply=blas \
  --disable-assert \
  CC=icc \
  CPPFLAGS="-I${MKLROOT}/include" \
  LDFLAGS="-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm"
