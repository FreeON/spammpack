#!/bin/bash

NUM_THREADS=( 48 24 12 1 2 4 8 16 20 28 32 36 40 44 )
#NUM_THREADS=( 40 44 )
TOL=( 1.0e-10 1.0e-8 1.0e-6 )

MATRIX=h2o_30.OrthoD
#MATRIX=h2o_90.OrthoD

for tol in ${TOL[@]}; do
  echo "tolerance = ${TOL}"
  for t in ${NUM_THREADS[@]}; do
    echo -n "$t numpy.mean([ "
    for (( i = 0; i < 5; i++ )); do
      OMP_NUM_THREADS=${t} \
        numactl --membind=0 \
        ./chunk_multiply \
        -T BCSR \
        -f ${MATRIX} \
        --tolerance ${TOL} \
        --no-verify \
        | grep seconds \
        | sed -e 's:^.*using \([0-9]\+\) OpenMP.* \([^ ]\+\) seconds:\1 \2:' \
        | awk '{printf("%s", $2);}'
      if [[ ${i} -lt 4 ]]; then
        echo -n ", "
      else
        echo " ])"
      fi
    done
  done
done
