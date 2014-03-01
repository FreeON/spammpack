#!/bin/bash

#./spamm-charm \
#  --block 256 \
#  --basic 4 \
#  --tolerance 0 \
#  --operation SP2 \
#  --fockian ~/data-sets/water/h2o_10.F_DIIS \
#  --density ~/data-sets/water/h2o_10.OrthoD \
#  --Ne 100 \
#  --iterations 100

./spamm-charm \
  --block 750 \
  --basic 750 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_30.F_DIIS \
  --density ~/data-sets/water/h2o_30.OrthoD \
  --Ne 300 \
  --iterations 100

./spamm-charm \
  --block 1250 \
  --basic 1250 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_50.F_DIIS \
  --density ~/data-sets/water/h2o_50.OrthoD \
  --Ne 500 \
  --iterations 100

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
