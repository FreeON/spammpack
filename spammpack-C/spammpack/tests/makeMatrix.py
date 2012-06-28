#!/usr/bin/env python

import argparse
import math
import random
import sys

parser = argparse.ArgumentParser()

parser.add_argument("-N",
    help = "Matrix dimension (default = %(default)s)",
    type = int,
    default = 1)

parser.add_argument("--type",
    dest = "matrixtype",
    help = "Matrix type (default = %(default)s)",
    choices = [ "full", "decay" ],
    default = "full")

parser.add_argument("--decay",
    metavar = "gamma",
    help = "Set matrix element decay to exp(-gamma |i-j|)",
    type = float,
    default = 1)

options = parser.parse_args()

A = [ [ 0 for i in range(options.N) ] for j in range(options.N) ]

if options.matrixtype == "full":
  for i in range(options.N):
    for j in range(options.N):
      A[i][j] = random.uniform(0.1, 1.0)

elif options.matrixtype == "decay":
  for i in range(options.N):
    A[i][i] = random.uniform(0.6, 1.0)
  for i in range(options.N):
    for j in range(options.N):
      if i != j:
        A[i][j] = A[i][i]*math.exp(-options.decay*math.fabs(i-j))

else:
  sys.stderr.write("[FIXME] unknown matrix type (%s)\n" % (options.matrixtype))
  sys.exit(1)

# Analyze matrix.
element_range = [ 0, 1e-8, 1e-6, 1e-4, 1, 10]

range_count = [ 0 for i in range(len(element_range)) ]
count_total = 0
min_element = 10
max_element = -10

for i in range(options.N):
  for j in range(options.N):
    if A[i][j] < min_element:
      min_element = A[i][j]
    if A[i][j] > max_element:
      max_element = A[i][j]

    for k in range(len(element_range)-1):
      if element_range[k] <= math.fabs(A[i][j]) and math.fabs(A[i][j]) < element_range[k+1]:
        range_count[k] += 1
        count_total += 1

for k in range(len(element_range)-1):
  sys.stderr.write("[%e,%e) = %d (%1.2f%%)\n" % (element_range[k],
      element_range[k+1], range_count[k],
      range_count[k]/float(count_total)*100))
sys.stderr.write("total = %d (%d)\n" % (count_total, options.N**2))
sys.stderr.write("min(A) = %e\n" % (min_element))
sys.stderr.write("max(A) = %e\n" % (max_element))

for i in range(options.N):
  for j in range(options.N):
    print("%d %d % 1.14e" % (i+1, j+1, A[i][j]))
