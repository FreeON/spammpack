import re
import numpy as np

def read_MM (filename = "", from_stdin = False):
  """Read a matrix in Matrix Market format and return a matrix object.
  """

  if from_stdin:
    fd = sys.stdin
  else:
    fd = open(filename)

  header = fd.readline()

  if re.search("^%%MatrixMarket\s+matrix\s+coordinate\s+double\s+general$", header) == None:
    raise Exception("incorrect file format")

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
