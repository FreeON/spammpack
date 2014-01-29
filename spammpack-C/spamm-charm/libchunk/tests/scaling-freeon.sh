#!/bin/bash

OPTIONS=(
    "-N 4096"
    "-b 64"
    "-T exp_decay"
    "-l 0.995"
    "-v"
    )

TOLERANCE=( 0 1e-8 1e-7 1e-6 )
THREADS=( 1 2 4 8 12 16 32 40 )

for tolerance in ${TOLERANCE[@]}; do
  # Get product complexity.
  numactl --interleave=all \
    ./chunk_multiply ${OPTIONS[*]} --tolerance ${tolerance} -c | tee --append scaling.output

  for P in ${THREADS[@]}; do
    OMP_NUM_THREADS=${P} \
      numactl --interleave=all -- \
      ./chunk_multiply ${OPTIONS[*]} --tolerance ${tolerance} \
      | tee --append scaling.output
  done
done

grep "done multiplying" scaling.output \
  | sed -e "s/^.* \([0-9]\+\) OpenMP.*, \(.*\) seconds/\1:\2,/"
