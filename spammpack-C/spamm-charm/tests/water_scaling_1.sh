#!/bin/bash

run_spamm() {
  local N=$1
  local N_basic=$2
  local tolerance=$3
  local Ne=$4
  local filename=$5

  ../src/spamm-charm \
    --block ${N} \
    --basic ${N_basic} \
    --tolerance ${tolerance} \
    --operation SP2 \
    --symbolic \
    --fockian ~/data-sets/water/${filename}.F_DIIS \
    --density ~/data-sets/water/${filename}.OrthoD \
    --Ne ${Ne} \
    --iterations 100
}

export OMP_NUM_THREADS=1

run_spamm  256 4 0      100 h2o_10
run_spamm  256 4 1.0e-9 100 h2o_10
run_spamm  256 4 1.0e-8 100 h2o_10
run_spamm  256 4 1.0e-7 100 h2o_10
run_spamm  256 4 1.0e-6 100 h2o_10
run_spamm  256 4 1.0e-5 100 h2o_10

run_spamm 1024 4 0      300 h2o_30
run_spamm 1024 4 1.0e-9 300 h2o_30
run_spamm 1024 4 1.0e-8 300 h2o_30
run_spamm 1024 4 1.0e-7 300 h2o_30
run_spamm 1024 4 1.0e-6 300 h2o_30
run_spamm 1024 4 1.0e-5 300 h2o_30

run_spamm 2048 4 0      500 h2o_50
run_spamm 2048 4 1.0e-9 500 h2o_50
run_spamm 2048 4 1.0e-8 500 h2o_50
run_spamm 2048 4 1.0e-7 500 h2o_50
run_spamm 2048 4 1.0e-6 500 h2o_50
run_spamm 2048 4 1.0e-5 500 h2o_50

run_spamm 2048 4 0      700 h2o_70
run_spamm 2048 4 1.0e-9 700 h2o_70
run_spamm 2048 4 1.0e-8 700 h2o_70
run_spamm 2048 4 1.0e-7 700 h2o_70
run_spamm 2048 4 1.0e-6 700 h2o_70
run_spamm 2048 4 1.0e-5 700 h2o_70

run_spamm 4096 4 0      900 h2o_90
run_spamm 4096 4 1.0e-9 900 h2o_90
run_spamm 4096 4 1.0e-8 900 h2o_90
run_spamm 4096 4 1.0e-7 900 h2o_90
run_spamm 4096 4 1.0e-6 900 h2o_90
run_spamm 4096 4 1.0e-5 900 h2o_90

run_spamm 4096 4 0      1100 h2o_110
run_spamm 4096 4 1.0e-9 1100 h2o_110
run_spamm 4096 4 1.0e-8 1100 h2o_110
run_spamm 4096 4 1.0e-7 1100 h2o_110
run_spamm 4096 4 1.0e-6 1100 h2o_110
run_spamm 4096 4 1.0e-5 1100 h2o_110

run_spamm 4096 4 0      1300 h2o_130
run_spamm 4096 4 1.0e-9 1300 h2o_130
run_spamm 4096 4 1.0e-8 1300 h2o_130
run_spamm 4096 4 1.0e-7 1300 h2o_130
run_spamm 4096 4 1.0e-6 1300 h2o_130
run_spamm 4096 4 1.0e-5 1300 h2o_130

run_spamm 4096 4 0      1500 h2o_150
run_spamm 4096 4 1.0e-9 1500 h2o_150
run_spamm 4096 4 1.0e-8 1500 h2o_150
run_spamm 4096 4 1.0e-7 1500 h2o_150
run_spamm 4096 4 1.0e-6 1500 h2o_150
run_spamm 4096 4 1.0e-5 1500 h2o_150
