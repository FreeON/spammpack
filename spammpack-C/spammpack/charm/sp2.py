#!/usr/bin/env python

## @file
#
# Run the SP2 algorithm on a Fockian matrix.
#
# @author Nicolas Bock <nicolas.bock@freeon.org>
# @author Matt Challacombe <matt.challacombe@freeon.org>

import argparse
import math
import numpy as np
import re
import sys

from spectral_bounds import *
from matrix_market import *

## Calculate the product complexity.
#
# @param P The matrix.
# @param tolerance The SpAMM tolerance.
# @param blocksize The hypothetical blocksize.
#
# @return The complexity.
def get_complexity (P, tolerance, blocksize):
  if blocksize <= 0:
    return 0

  N = 1
  while N*blocksize < P.shape[0]:
    N *= 2
  #print("{:d}x{:d} matrix -> N = {:d}".format(P.shape[0], P.shape[1], N))

  norm = np.matrix(np.zeros(( N, N )))
  for i in range(N):
    for j in range(N):
      norm[i, j] = np.linalg.norm(P[(i*N):((i+1)*N),(j*N):((j+1)*N)])

  complexity = 0
  for i in range(N):
    for j in range(N):
      for k in range(N):
        if norm[i, k]*norm[k, j] > tolerance:
          complexity += 1

  return complexity, N**3

## Print a histogram of the matrix elements.
#
# @param P The matrix.
def print_histogram (P):
  print("matrix element magnitudes:")
  hist, bin_edges = np.histogram(abs(P[abs(P) > 0]),
      [ 0, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2,
        1e-1, 1.0, 2.0 ])
  count_width = int(math.ceil(math.log(P.shape[0]*P.shape[1], 10)))
  format_string = "[ {:1.2e}, {:1.2e} ) = {:" + str(count_width) + "d} {:" + str(count_width) + "d} {:6.3f}%"
  cumulative = 0
  for i in range(len(bin_edges)-1):
    cumulative += hist[i]
    print(format_string.format(
      bin_edges[i], bin_edges[i+1], hist[i], cumulative,
      100*cumulative/P.shape[0]/P.shape[1]))

## The main program.
def main ():
  parser = argparse.ArgumentParser()

  parser.add_argument("FILE",
      help = "The Fockian matrix, in matrix market format.",
      nargs = "*")

  parser.add_argument("--Ne",
      help = "The number of electrons (default %(default)s)",
      type = int,
      default = 1)

  parser.add_argument("--max-iterations",
      help = "The maximum number of iterations (default %(default)s)",
      metavar = "MAX",
      type = int,
      default = 100)

  parser.add_argument("--bin-density",
      help = "Print histogram of binned matrix elements of density guess",
      action = "store_true",
      default = False)

  parser.add_argument("--tolerance",
      help = "The hypothetical SpAMM tolerance",
      type = float,
      default = 0)

  parser.add_argument("--blocksize",
      help = "Set a hypothetical blocksize to calculate product complexities",
      type = int,
      default = 0)

  parser.add_argument("--debug", "-d",
      help = "Print debugging stuff",
      action = "store_true",
      default = False)

  options = parser.parse_args()

  F = None
  if len(options.FILE) == 0:
    F = read_MM()
  else:
    F = read_MM(options.FILE[0])

  bounds = spectral_bounds(F)
  print("spectral bounds of F: " + str(bounds))

  P = (bounds[1]*np.matrix(np.identity(F.shape[0]))-F)/(bounds[1]-bounds[0])
  print("spectral bounds of P: " + str(spectral_bounds(P)))

  if options.debug:
    print(P)
    write_MM(P)

  if options.bin_density:
    print_histogram(P)

  if options.max_iterations < 20:
    i_min = options.max_iterations
  else:
    i_min = 20

  occupation = np.zeros(4)
  trace_P = np.trace(P)
  print("iteration  0: trace(P) = {:e}".format(trace_P))

  converged = False
  for i in range(options.max_iterations):
    P2 = P*P
    trace_P2 = np.trace(P2)
    complexity, full_complexity = get_complexity(P, options.tolerance, options.blocksize)

    if options.debug:
      print(P2)

    if abs(trace_P2-options.Ne/2.0) < abs(2*trace_P-trace_P2-options.Ne/2.0):
      P = P2
      trace_P = trace_P2
    else:
      P = 2*P-P2
      trace_P = 2*trace_P-trace_P2

    occupation[0] = trace_P

    print("iteration {:2d}: ".format(i+1)
        + "trace(P) = {:e}, ".format(trace_P)
        + "trace(P)-Ne/2 = {: e}, ".format(trace_P-options.Ne/2.0)
        + "complexity = {:d} ".format(complexity)
        + "(out of {:d})".format(full_complexity))

    if i > i_min:
      current_error = abs(occupation[1]-occupation[0])
      if current_error < 1e-2:
        last_error = abs(occupation[3]-occupation[2])
        if last_error <= current_error:
          converged = True
          break

    for i in range(3, 0, -1):
      occupation[i] = occupation[i-1]

  if converged:
    print("converged in {:d} steps".format(i+1))

    if options.bin_density:
      print_histogram(P)

  else:
    print("failed to converge")
    sys.exit(1)

##############################################

if __name__ == "__main__":
  main()
