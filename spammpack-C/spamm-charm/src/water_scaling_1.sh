#!/bin/bash

run_spamm() {
  local N=$1
  local N_basic=$2
  local tolerance=$3
  local Ne=$4
  local filename=$5

  ./spamm-charm \
    --block ${N} \
    --basic ${N_basic} \
    --tolerance ${tolerance} \
    --operation SP2 \
    --fockian ~/data-sets/water/${filename}.F_DIIS \
    --density ~/data-sets/water/${filename}.OrthoD \
    --Ne ${Ne} \
    --iterations 100
}

export OMP_NUM_THREADS=1

run_spamm  256 4 0      100 h2o_10
run_spamm  256 4 1.0e-8 100 h2o_10
run_spamm  256 4 1.0e-7 100 h2o_10
run_spamm  256 4 1.0e-6 100 h2o_10
run_spamm  256 4 1.0e-5 100 h2o_10

run_spamm 1024 4 0      300 h2o_30
run_spamm 1024 4 1.0e-8 300 h2o_30
run_spamm 1024 4 1.0e-7 300 h2o_30
run_spamm 1024 4 1.0e-6 300 h2o_30
run_spamm 1024 4 1.0e-5 300 h2o_30

run_spamm 2048 4 0      500 h2o_50
run_spamm 2048 4 1.0e-8 500 h2o_50
run_spamm 2048 4 1.0e-7 500 h2o_50
run_spamm 2048 4 1.0e-6 500 h2o_50
run_spamm 2048 4 1.0e-5 500 h2o_50

exit 0

./spamm-charm \
  --block 1750 \
  --basic 1750 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_70.F_DIIS \
  --density ~/data-sets/water/h2o_70.OrthoD \
  --Ne 700 \
  --iterations 100

./spamm-charm \
  --block 2250 \
  --basic 2250 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_90.F_DIIS \
  --density ~/data-sets/water/h2o_90.OrthoD \
  --Ne 900 \
  --iterations 100

./spamm-charm \
  --block 2750 \
  --basic 2750 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_110.F_DIIS \
  --density ~/data-sets/water/h2o_110.OrthoD \
  --Ne 1100 \
  --iterations 100

./spamm-charm \
  --block 3250 \
  --basic 3250 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_130.F_DIIS \
  --density ~/data-sets/water/h2o_130.OrthoD \
  --Ne 1300 \
  --iterations 100

./spamm-charm \
  --block 3750 \
  --basic 3750 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_150.F_DIIS \
  --density ~/data-sets/water/h2o_150.OrthoD \
  --Ne 1500 \
  --iterations 100
