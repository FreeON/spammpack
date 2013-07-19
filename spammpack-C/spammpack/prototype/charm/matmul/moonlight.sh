#!/bin/bash

srun --mail-type=ALL --nodes=1 ./charmrun -np 1 matmul -N 2048 -b 64
