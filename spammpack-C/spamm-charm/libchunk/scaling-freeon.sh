#!/bin/bash

#./configure-mustang.sh parallel
#make clean && make
#bash -c "cd tests && ./scaling-freeon-dense.sh"

./configure-mustang.sh sequential
make clean && make
bash -c "cd tests && ./scaling-freeon-SpAMM.sh"
