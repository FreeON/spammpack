#!/bin/bash

for i in 20; do
  sed -e "s:nodes=[0-9]\+:nodes=${i}:" \
    -e "s:np\s\+[0-9]\+:np ${i}:" \
    -e "s:matmul-[0-9]\+:matmul-${i}:" \
    moonlight.sh.template > moonlight.sh
  sbatch moonlight.sh
done
