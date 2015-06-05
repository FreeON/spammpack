#!/bin/bash

./configure CC=gcc CFLAGS="-Ofast -g" --disable-shared --disable-reset-C --enable-chunk-tree-max-tier=8
make clean
make

LOGFILE=scan-chunk-sizes

for N_chunk in 2048 1024 512 256 128; do
  for N_basic in 128 64 32 16; do
    echo "# N_chunk = ${N_chunk}, N_basic = ${N_basic}" >> ${LOGFILE}.dat
    for P in 1 16 24 48; do
      export OMP_NUM_THREADS=${P}
      ./tests/chunk_multiply \
        --type full \
        -N 2048 \
        --N_chunk ${N_chunk} \
        --N_basic ${N_basic} \
        --repeat 5 \
        --no-verify
    done \
      | tee --append ${LOGFILE}.output \
      | awk '/products/ { if(!header) { printf("# %d products\n", $9); header = 1; } } /seconds/ { if(!n) { n = 0; }; P[n] = $1; time[n] = $15; std[n] = $17; printf("% 3d %e %e %6.2f %6.2f%%\n", P[n], time[n], std[n], time[0]/time[n], time[0]/time[n]/P[n]*100); n++; }' \
      | tee --append ${LOGFILE}.dat
  done
done
