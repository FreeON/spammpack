#!/bin/bash

. ./scaling-freeon-common.sh

TOLERANCE=( 0 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8 )

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
