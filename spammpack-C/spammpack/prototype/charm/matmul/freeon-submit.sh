#!/bin/bash

for i in 1 2 3 8 12 16 20; do
  sed \
    -e "s:procs=[0-9]\+:procs=${i}:" \
    -e "s:NUMBER_NODES=[0-9]\+:NUMBER_NODES=${i}:" \
    freeon.job.template > freeon.job
  qsub freeon.job
done
