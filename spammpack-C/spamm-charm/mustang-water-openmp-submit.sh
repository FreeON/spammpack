#!/bin/bash

run_job() {
  local NAME=$1
  local NODES=$2
  local NE=$3
  local TOLERANCE=$4
  local SPAMM_VERSION=$5
  local BLOCK=$6
  local BASIC=$7

  sed \
  -e "s:nodes=N:nodes=${NODES}:" \
  -e "s:JOBNAME:${NAME}-${SPAMM_VERSION}:" \
  -e "s:SPAMM_VERSION:${SPAMM_VERSION}:" \
  -e "s:NE:${NE}:" \
  -e "s:FOCKIAN:/scratch/nbock/${NAME}.F_DIIS:" \
  -e "s:TOLERANCE:${TOLERANCE}:" \
  -e "s:BLOCK:${BLOCK}:" \
  -e "s:BASIC:${BASIC}:" \
  mustang.job.water.template | tee mustang.job
  msub mustang.job
}

#NODES=( 1 4 12 32 64 128 256 )
NODES=( 768 1024 )
WATERS=( h2o_10 h2o_30 h2o_50 h2o_70 h2o_90 h2o_110 h2o_130 h2o_150 )
NE=( 100 300 500 700 900 1100 1300 1500 )
TOLERANCES=( 1e-10 1e-8 1e-6 )
BLOCK=128
BASIC=4
SPAMM_VERSION=openmp

run_job h2o_50 1 500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_50 1 500 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_50 1 500 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_50 4 500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}

run_job h2o_70  1 700 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_70  1 700 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_70  1 700 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_70  4 700 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_70  4 700 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_70  4 700 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_70 12 700 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}

run_job h2o_90  1 900 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90  1 900 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90  1 900 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90  4 900 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90  4 900 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90  4 900 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90 12 900 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_90 12 900 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}

run_job h2o_110  1 1100 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110  1 1100 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110  1 1100 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110  4 1100 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110  4 1100 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110  4 1100 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110 12 1100 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110 12 1100 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_110 32 1100 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}

run_job h2o_130  1 1300 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130  1 1300 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130  1 1300 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130  4 1300 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130  4 1300 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130  4 1300 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130 12 1300 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130 12 1300 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130 12 1300 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130 32 1300 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_130 32 1300 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}

run_job h2o_150  1 1500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150  1 1500 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150  1 1500 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150  4 1500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150  4 1500 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150  4 1500 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150 12 1500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150 12 1500 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150 12 1500 1e-06 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150 32 1500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150 32 1500 1e-08 ${SPAMM_VERSION} ${BLOCK} ${BASIC}
run_job h2o_150 64 1500 1e-10 ${SPAMM_VERSION} ${BLOCK} ${BASIC}

exit

for i in ${NODES[@]}; do
  for j in ${TOLERANCES[@]}; do
    for (( k = 0; k < ${#WATERS[@]}; k++ )); do
      run_job ${WATERS[$k]} ${i} ${NE[$k]} ${j} ${SPAMM_VERSION} ${BLOCK} ${BASIC}
    done
  done
done
rm -f mustang.job
