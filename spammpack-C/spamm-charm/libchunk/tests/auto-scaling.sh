#!/bin/bash

for BUILD_TYPE in serial openmp; do
  for BUILD_COMPILER in gcc intel; do
    ./scaling-freeon.sh
  done
done
