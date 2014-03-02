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
    return 0, 0

  N = 1
  while N*blocksize < P.shape[0]:
    N *= 2

  PPadded = np.matrix(np.zeros(( N*blocksize, N*blocksize )))
  PPadded[0:P.shape[0], 0:P.shape[1]] = P

  norm = np.matrix(np.zeros(( N, N )))
  for i in range(N):
    for j in range(N):
      norm[i, j] = np.linalg.norm(PPadded[(i*blocksize):((i+1)*blocksize),
        (j*blocksize):((j+1)*blocksize)], "fro")

  complexity = 0
  for i in range(N):
    for j in range(N):
      for k in range(N):
        if norm[i, k]*norm[k, j] > tolerance:
          #print("keeping ({:d},{:d},{:d}) ".format(i, j, k)
          #    + "A[{:d}:{:d},{:d}:{:d}]*".format(
          #      i*blocksize, (i+1)*blocksize, k*blocksize, (k+1)*blocksize)
          #    + "B[{:d}:{:d},{:d}:{:d}] ".format(
          #      k*blocksize, (k+1)*blocksize, j*blocksize, (j+1)*blocksize)
          #    + "({:e} * {:e} = {:e} ".format(
          #      norm[i, k], norm[k, j], norm[i, k]*norm[k, j]))
          complexity += 1
        #else:
        #  print("pruning ({:d},{:d},{:d}) ".format(i, j, k)
        #      + "A[{:d}:{:d},{:d}:{:d}]*".format(
        #        i*blocksize, (i+1)*blocksize, k*blocksize, (k+1)*blocksize)
        #      + "B[{:d}:{:d},{:d}:{:d}] ".format(
        #        k*blocksize, (k+1)*blocksize, j*blocksize, (j+1)*blocksize)
        #      + "({:e} * {:e} = {:e} ".format(
        #        norm[i, k], norm[k, j], norm[i, k]*norm[k, j]))

  return complexity, (math.ceil(P.shape[0]/blocksize))**3

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

  parser.add_argument(
      "FILE",
      help = "The Fockian matrix, in matrix market format."
      )

  parser.add_argument(
      "--density",
      help = "The converged density for verification "
      + "(in matrix market format)."
      )

  parser.add_argument(
      "--Ne",
      help = "The number of electrons (default %(default)s)",
      type = int,
      default = 1
      )

  parser.add_argument(
      "--max-iterations",
      help = "The maximum number of iterations (default %(default)s)",
      metavar = "MAX",
      type = int,
      default = 100
      )

  parser.add_argument(
      "--bin-density",
      help = "Print histogram of binned matrix elements of density guess",
      action = "store_true",
      default = False
      )

  parser.add_argument(
      "--tolerance",
      help = "The hypothetical SpAMM tolerance",
      type = float,
      default = 0
      )

  parser.add_argument(
      "--blocksize",
      help = "Set a hypothetical blocksize to calculate product complexities",
      type = int,
      default = 0
      )

  parser.add_argument(
      "--debug", "-d",
      help = "Print debugging stuff",
      action = "store_true",
      default = False
      )

  options = parser.parse_args()

  F = read_MM(options.FILE)

  if options.debug:
    for i in range(F.shape[0]):
      for j in range(F.shape[1]):
        print("{:d} {:d} {: e}".format(i, j, F[i, j]))

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
  print("iteration  0: trace(P) = {:e}, ||P|| = {:1.16e}".format(
    trace_P,
    np.linalg.norm(P, "fro")))

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
        + "||P|| = {:1.16e}, ".format(np.linalg.norm(P, "fro"))
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
    print("SP2 converged in {:d} steps".format(i+1))
    print("idempotency error          = {:e}".format(abs(occupation[1]-occupation[0])))
    print("previous idempotency error = {:e}".format(abs(occupation[3]-occupation[2])))

    if options.bin_density:
      print_histogram(P)

    if options.density:
      print("loading reference density")
      D = read_MM(options.density)
      print("||D|| = {:1.16e}, ||P-D|| = {:1.16e}".format(
        np.linalg.norm(D, "fro"),
        np.linalg.norm(P-D, "fro")
        ))

  else:
    print("failed to converge")
    sys.exit(1)

##############################################

if __name__ == "__main__":
  main()
