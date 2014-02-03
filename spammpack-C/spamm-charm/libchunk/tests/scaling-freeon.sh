#!/bin/bash

#N=8192
#B=128
N=4096
B=32
REPEAT=10

OPTIONS=(
    "-N ${N}"
    "-b ${B}"
    "-T exp_decay"
    "-l 0.994"
    "-R ${REPEAT}"
    "-c"
    "-v"
    )

TOLERANCE=( 0 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 )
THREADS=( 48 1 2 4 8 12 16 20 24 28 32 36 40 )

for tolerance in ${TOLERANCE[@]}; do
  for P in ${THREADS[@]}; do
    OMP_NUM_THREADS=${P} \
      numactl --interleave=all -- \
      ./chunk_multiply ${OPTIONS[*]} --tolerance ${tolerance} \
      | tee --append scaling_OpenMP.output
  done
done

grep "done multiplying" scaling_OpenMP.output \
  | sed -e "s/^.* \([0-9]\+\) OpenMP.*, \(.*\) seconds/\1:\2,/"
