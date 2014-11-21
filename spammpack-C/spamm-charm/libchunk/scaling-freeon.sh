#!/bin/bash

./configure-freeon.sh parallel
make clean && make
bash -c "cd tests && ./scaling-freeon-dense.sh"

./configure-freeon.sh sequential
make clean && make
bash -c "cd tests && ./scaling-freeon-SpAMM.sh"
