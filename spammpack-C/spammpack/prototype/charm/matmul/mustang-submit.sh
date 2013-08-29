#!/bin/bash

NODES=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 96 128 160 192 224 256 )

for i in ${NODES[@]}; do

  sed -e "s:nodes=[0-9]\+:nodes=${i}:" \
    -e "s:MM-[0-9]\+:MM-${i}:" \
    mustang.job.full.template > mustang.job
  msub mustang.job

  sed -e "s:nodes=[0-9]\+:nodes=${i}:" \
    -e "s:MM-[0-9]\+:MM-${i}:" \
    mustang.job.decay.template > mustang.job
  msub mustang.job

  sed -e "s:nodes=[0-9]\+:nodes=${i}:" \
    -e "s:MM-[0-9]\+:MM-${i}:" \
    mustang.job.diagonal.template > mustang.job
  msub mustang.job

done
rm -f mustang.job
