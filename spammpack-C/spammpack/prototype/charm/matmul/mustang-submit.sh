#!/bin/bash

NODES=( 1 2 4 8 12 16 20 24 28 32 36 40 44 48 52 56 60 )

for i in ${NODES[@]}; do
  sed -e "s:nodes=[0-9]\+:nodes=${i}:" \
    -e "s:MM-[0-9]\+:MM-${i}:" \
    mustang.job.template > mustang.job
  msub mustang.job
  rm -f mustang.job
done
