#!/usr/bin/python

import random, optparse

parser = optparse.OptionParser()

parser.add_option("-N",
    help = "generate NxN matrix",
    type = "int",
    default = 2,
    dest = "N",
    metavar = "N")

(options, arguments) = parser.parse_args()

print "%%MatrixMarket matrix coordinate real general"
print "%i %i %i" % (options.N, options.N, options.N*options.N)
for i in range(options.N):
  for j in range(options.N):
    print "%i %i %f" % (i+1, j+1, random.random()+2)
