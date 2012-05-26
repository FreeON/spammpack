#!/usr/bin/env python

import argparse
import random

parser = argparse.ArgumentParser()

parser.add_argument("-N",
    help = "Matrix dimension (default = %(default)s)",
    type = int,
    default = 1)

options = parser.parse_args()

for i in range(options.N):
  for j in range(options.N):
    print "%d %d % 1.14e" % (i+1, j+1, random.uniform(0.1, 1.0))
