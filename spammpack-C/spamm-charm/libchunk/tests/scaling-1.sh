#!/bin/bash
#
# Building with gcc:
pushd ..
./configure --disable-assert --enable-chunk-tree-max-tier=4 CFLAGS=-"O3 -ftree-vectorize -g"
make clean
make
popd

for tolerance in 1e-10 1e-6; do
  for threads in 1 2 4 8 12 16 20 24 28 32 36 40 44 48; do
    OMP_NUM_THREADS=${threads} numactl \
      --interleave=all -- \
      ./chunk_multiply \
      --tolerance ${tolerance} \
      --N_basic 4 \
      --type BCSR \
      --bcsr ~/data-sets/water/h2o_90.OrthoD \
      --no-verify \
      --repeat 10 | tee --append OpenMP.output
  done
done
