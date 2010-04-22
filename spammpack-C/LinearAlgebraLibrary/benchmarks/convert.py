#!/usr/bin/python

import math, sys

fd = open("SparseHamiltonianNick.dat")
lines = fd.readlines()
fd.close()

print >> sys.stderr, "read input, converting to MM format"

number_dropped = 0
line_number = 0
min_element = 100
max_element = 0
for i in range(1600):
  for j in range(1600):
    x = float(lines[line_number])
    line_number += 1
    if math.fabs(x) > 1e-14:
      print i+1, j+1, x
      if math.fabs(x) < min_element:
        min_element = math.fabs(x)
      if math.fabs(x) > max_element:
        max_element = math.fabs(x)
    else:
      number_dropped += 1

print >> sys.stderr, number_dropped
print >> sys.stderr, "max element = ", max_element
print >> sys.stderr, "min element = ", min_element
