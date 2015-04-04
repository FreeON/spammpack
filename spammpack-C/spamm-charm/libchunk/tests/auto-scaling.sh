#!/bin/bash

for BUILD_TYPE in serial openmp; do
  for BUILD_COMPILER in gcc intel; do
    BUILD_TYPE=${BUILD_TYPE} BUILD_COMPILER=${BUILD_COMPILER} ./scaling-freeon.sh
  done
done
