#!/bin/bash

./spamm-charm \
  --block 512 \
  --basic 4 \
  --tolerance 0 \
  --operation SP2 \
  --fockian ~/data-sets/water/h2o_10.OrthoF \
  --density ~/data-sets/water/h2o_10.OrthoD \
  --Ne 100 \
  --F-min -1.919805e+01 \
  --F-max 4.192886e+00
