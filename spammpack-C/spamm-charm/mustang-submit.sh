#!/bin/bash

NODES=( 1 2 4 16 32 64 128 256 )

DENSITY=( PP020-2 )
NE=( 1082 )
BLOCK=( 256 )
BASIC=( 16 )

for i in ${NODES[@]}; do
  for (( j = 0; j < ${#DENSITY[@]}; j++ )); do

    sed \
      -e "s:nodes=N:nodes=${i}:" \
      -e "s:JOBNAME:${DENSITY[$j]}:" \
      -e "s:NE:${NE[$j]}:" \
      -e "s:DENSITY:~/Fockians/${DENSITY[$j]}.OrthoF:" \
      -e "s:PEMAP:${DENSITY[$j]}-${i}:" \
      -e "s:BLOCK:${BLOCK[$j]}:" \
      -e "s:BASIC:${BASIC[$j]}:" \
      mustang.job.PP.template > mustang.job
    msub mustang.job

  done
done
rm -f mustang.job
