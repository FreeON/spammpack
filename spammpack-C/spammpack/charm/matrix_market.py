## @file
#
# Some matrix market functions.
#
# @author Nicolas Bock <nicolas.bock@freeon.org>
# @author Matt Challacombe <matt.challacombe@freeon.org>

import re
import sys
import numpy as np

## Read a matrix in Matrix Market format and return a matrix object.
def read_MM (filename = ""):
  if len(filename) == 0:
    print("reading from stdin")
    fd = sys.stdin
  else:
    print("reading from \"" + filename + "\"")
    fd = open(filename)

  header = fd.readline()

  if re.search("^%%MatrixMarket\s+matrix\s+coordinate\s+double\s+general$", header) == None:
    raise Exception("incorrect file header")

  comment = re.compile("^%")

  while True:
    line = fd.readline()
    if comment.search(line):
      continue
    else:
      break

  token = line.split()
  M = int(token[0])
  N = int(token[1])
  nnz = int(token[2])

  F = np.matrix(np.zeros((M, N), dtype = 'float'))

  for line in fd:
    token = line.split()
    F[int(token[0])-1, int(token[1])-1] = float(token[2])

  fd.close()

  return F

## Write a matrix in Matrix Market format to standard output if the filename
# argument is missing or to a file.
def write_MM (A, filename = ""):
  print("%%MatrixMarket matrix coordinate double general")
  print(A.shape[0], A.shape[1], A.shape[0]*A.shape[1])
  for i in range(A.shape[0]):
    for j in range(A.shape[1]):
      print("{:d} {:d} {: e}".format(i+1, j+1, A[i,j]))
