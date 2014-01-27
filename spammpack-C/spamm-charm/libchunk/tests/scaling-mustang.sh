#!/bin/bash

for P in 24 1 2 4 8 12 16 20; do
  OMP_NUM_THREADS=$P \
    numactl --interleave=all -- \
    ./chunk_multiply -N 4096 -b 64 -v \
    | tee --append scaling.output
done

grep "done multiplying" scaling.output \
  | sed -e "s/^.* \([0-9]\+\) OpenMP.*, \(.*\) seconds/\1:\2,/"
