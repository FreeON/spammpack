#!/usr/bin/env python

import argparse
import numpy
import re
import sys

##############################################

def add2D (i, j, PE, PEMap):
  if not i in PEMap:
    PEMap[i] = {}
  PEMap[i][j] = PE

def printMap (PEMap):
  currentI = -1
  currentJ = -1

##############################################

parser = argparse.ArgumentParser()

parser.add_argument("FILE",
    help = "The output file to plot. A value of '-' means standard input.")

parser.add_argument("--debug",
    help = "Print debugging stuff",
    action = "store_true",
    default = False)

options = parser.parse_args()

if options.FILE == "-":
  fd = sys.stdin
else:
  fd = open(options.FILE)

inIteration = False
iteration = -1

inMap = False
currentMap = ""

PEMap = {}

for line in fd:
  if options.debug:
    print("read: ", line.rstrip())

  if not inIteration:
    result = re.compile("iteration ([0-9]+) on").search(line)
    if result:
      inIteration = True
      iteration = int(result.group(1))
      print("iteration {:d}".format(iteration))
      continue

  result = re.compile("] PEMap for (.*):").search(line)
  if result:
    mapName = result.group(1)
    if inMap and mapName != currentMap:
      raise(Exception("map {:s} already open for reading".format(mapName)))
    if not inMap:
      N = 0
      elementBuffer = []
    inMap = True
    currentMap = mapName
    if options.debug:
      print("opening map {:s}".format(currentMap))

  result = re.compile("PEMap\(([0-9]+),([0-9]+)\) = ([0-9]+)").search(line)
  if result:
    i = int(result.group(1))
    j = int(result.group(2))
    PE = int(result.group(3))
    if not inMap:
      raise(Exception("no map open for reading"))
    elementBuffer.append( (i, j, PE) )
    if i+1 > N:
      N = i+1
    if j+1 > N:
      N = j+1
    continue

  result = re.compile("PEMap\(([0-9]+),([0-9]+),([0-9]+)\) = ([0-9]+)").search(line)
  if result:
    i = int(result.group(1))
    j = int(result.group(2))
    k = int(result.group(3))
    PE = int(result.group(4))
    if not inMap:
      raise(Exception("no map open for reading"))
    elementBuffer.append( (i, j, k, PE) )
    if i+1 > N:
      N = i+1
    if j+1 > N:
      N = j+1
    if k+1 > N:
      N = k+1
    continue

  result = re.compile("end of PEMap").search(line)
  if result:
    if currentMap == "convolution":
      PEMap[currentMap] = numpy.empty([N, N, N])
      for (i, j, k, PE) in elementBuffer:
        PEMap[currentMap][i,j,k] = PE
    else:
      PEMap[currentMap] = numpy.empty([N, N])
      for (i, j, PE) in elementBuffer:
        PEMap[currentMap][i,j] = PE
    print(PEMap[currentMap])
    inMap = False
    if options.debug:
      print("closing map {:s}".format(currentMap))
    continue

fd.close()

for name in PEMap:
  printMap(PEMap[name])
