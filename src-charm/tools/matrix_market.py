## @file
#
# Some matrix market functions.
#
# @author Nicolas Bock <nicolas.bock@freeon.org>
# @author Matt Challacombe <matt.challacombe@freeon.org>

## Read a matrix in Matrix Market format and return a matrix object.
def read_MM (filename = ""):
  import numpy as np
  import re
  import sys

  if len(filename) == 0:
    sys.stdout.write("reading from stdin... ")
    fd = sys.stdin
  else:
    sys.stdout.write("reading from \"" + filename + "\"... ")
    fd = open(filename)
  sys.stdout.flush()

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

  print("||F|| = {:1.16e}".format(np.linalg.norm(F, "fro")))

  return F

## Write a matrix in Matrix Market format to standard output if the filename
# argument is missing or to a file.
def write_MM (A, filename = ""):
  print("%%MatrixMarket matrix coordinate double general")
  print(A.shape[0], A.shape[1], A.shape[0]*A.shape[1])
  for i in range(A.shape[0]):
    for j in range(A.shape[1]):
      print("{:d} {:d} {: e}".format(i+1, j+1, A[i,j]))
