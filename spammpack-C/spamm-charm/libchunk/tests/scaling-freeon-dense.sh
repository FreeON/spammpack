#!/bin/bash

. ./scaling-freeon-common.sh

for P in ${THREADS[@]}; do
  OMP_NUM_THREADS=${P} \
    numactl --interleave=all -- \
    ./chunk_multiply ${OPTIONS[*]} --dense \
    | tee --append scaling_OpenMP.output
done

grep "done multiplying" scaling_OpenMP.output \
  | sed -e "s/^.* \([0-9]\+\) OpenMP.*, \(.*\) seconds/\1:\2,/"
