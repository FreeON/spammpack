#!/bin/bash

for i in 1 2 4 8 12 16 20 24 28 32 36 40; do
  sed -e "s:nodes=[0-9]\+:nodes=${i}:" \
    -e "s:matmul-[0-9]\+:matmul-${i}:" \
    mustang.job.template > mustang.job
  msub mustang.job
  rm -f mustang.job
done
