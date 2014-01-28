#!/bin/bash

OPTIONS=(
    "-N 4096"
    "-b 64"
    "-T exp_decay"
    "-l 0.995"
    "-t 1.0e-8"
    "-v"
    )

# Get product complexity.
./chunk_multiply ${OPTIONS[*]} -c | tee --append scaling.output

for P in 1 2 4 8 12 16 32 20 24; do
  OMP_NUM_THREADS=$P \
    numactl --interleave=all -- \
    ./chunk_multiply ${OPTIONS[*]} \
    | tee --append scaling.output
done

grep "done multiplying" scaling.output \
  | sed -e "s/^.* \([0-9]\+\) OpenMP.*, \(.*\) seconds/\1:\2,/"
