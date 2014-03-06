#!/bin/bash

NODES=( 1 )
WATERS=( h2o_10 h2o_30 h2o_50 h2o_70 h2o_90 h2o_110 h2o_130 h2o_150 )
NE=( 100 300 500 700 900 1100 1300 1500 )
TOLERANCES=( 1e-10 1e-8 1e-6 )
BLOCK=128
BASIC=4
SPAMM_VERSION=serial

for i in ${NODES[@]}; do
  for j in ${TOLERANCES[@]}; do
    for (( k = 0; k < ${#WATERS[@]}; k++ )); do

      sed \
      -e "s:nodes=N:nodes=${i}:" \
      -e "s:JOBNAME:${WATERS[$k]}-${SPAMM_VERSION}:" \
      -e "s:SPAMM_VERSION:${SPAMM_VERSION}:" \
      -e "s:NE:${NE[$k]}:" \
      -e "s:FOCKIAN:/scratch/nbock/${WATERS[$k]}.F_DIIS:" \
      -e "s:TOLERANCE:${j}:" \
      -e "s:BLOCK:${BLOCK}:" \
      -e "s:BASIC:${BASIC}:" \
      mustang.job.water.template | tee mustang.job
      msub mustang.job

    done
  done
done
rm -f mustang.job
