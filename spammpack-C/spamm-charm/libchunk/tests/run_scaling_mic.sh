#!/bin/bash

if [[ $# -ne 2 ]]; then
  echo "missing host name and directory"
  exit 1
fi
MIC=$1
MIC_DIR=$2

NUM_THREADS=( 240 120 60 30 1 )
#NUM_THREADS=( 40 44 )
TOL=( 1.0e-10 1.0e-8 1.0e-6 )

#MATRIX=h2o_30.OrthoD
MATRIX=h2o_90.OrthoD

# Test with full output.
ssh $MIC "cd $MIC_DIR; LD_LIBRARY_PATH=${MIC_DIR} ./chunk_multiply --type BCSR --bcsr ${MATRIX} --tolerance 0.0e0 --no-verify"

for tol in ${TOL[@]}; do
  echo "tolerance = ${tol}"
  for t in ${NUM_THREADS[@]}; do
    echo -n "$t; numpy.mean([ "
    for (( i = 0; i < 5; i++ )); do
      ssh $MIC "cd $MIC_DIR; LD_LIBRARY_PATH=${MIC_DIR} OMP_NUM_THREADS=${t} ./chunk_multiply --type BCSR --bcsr ${MATRIX} --tolerance ${tol} --no-verify" \
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
