#!/bin/bash

./spamm-charm \
  --block 250 \
  --basic 250 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_10.F_DIIS \
  --density ~/data-sets/water/h2o_10.OrthoD \
  --Ne 100 \
  --iterations 100

#  --F-min -1.919805e+01
#  --F-max 4.192886e+00
