#!/usr/bin/env python

import argparse
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

options = parser.parse_args()

if options.matrixtype == "full":
  for i in range(options.N):
    for j in range(options.N):
      print "%d %d % 1.14e" % (i+1, j+1, random.uniform(0.1, 1.0))

elif options.matrixtype == "decay":
  print "not implemented"
  sys.exit(1)

else:
  print "unknown matrix type (%s)" % (options.matrixtype)
  sys.exit(1)
