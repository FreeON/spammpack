#!/bin/bash

NODES=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 96 128 160 192 224 256 )

DENSITY=( PP040 PP080 PP120 )
NE=( 2162 4322 6482 )
BLOCK=( 128 128 128 )

for i in ${NODES[@]}; do
  for (( j = 0; j < ${#DENSITY[@]}; j++ )); do

    sed \
      -e "s:nodes=N:nodes=${i}:" \
      -e "s:JOBNAME:${DENSITY[$j]}:" \
      -e "s:NE:${NE[$j]}:" \
      -e "s:DENSITY:~/Fockians/${DENSITY[$j]}.OrthoF:" \
      -e "s:PEMAP:${DENSITY[$j]}-${i}:" \
      -e "s:BLOCK:${BLOCK[$j]}:" \
      mustang.job.PP.template > mustang.job
    msub mustang.job

  done
done
rm -f mustang.job
