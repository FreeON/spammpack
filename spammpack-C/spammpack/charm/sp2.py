#!/usr/bin/env python

## @file
#
# Run the SP2 algorithm on a Fockian matrix.
#
# @author Nicolas Bock <nicolas.bock@freeon.org>
# @author Matt Challacombe <matt.challacombe@freeon.org>

import argparse
import numpy as np
import re
import sys

from spectral_bounds import *
from matrix_market import *

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
      help = "Print histogram of binned matrix elements of density guess\n",
      action = "store_true",
      default = False)

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
    print("matrix element magnitudes:")
    hist, bin_edges = np.histogram(abs(P[abs(P) > 0]),
        [ 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0 ])
    cumulative = 0
    for i in range(len(bin_edges)-1):
      cumulative += hist[i]
      print("[ {:e}, {:e} ) = {:5d} {:5d} {:6.3f}%".format(
        bin_edges[i], bin_edges[i+1], hist[i], cumulative,
        100*cumulative/P.shape[0]/P.shape[1]))

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

    if options.debug:
      print(P2)

    if abs(trace_P2-options.Ne/2.0) < abs(2*trace_P-trace_P2-options.Ne/2.0):
      P = P2
      trace_P = trace_P2
    else:
      P = 2*P-P2
      trace_P = 2*trace_P-trace_P2

    occupation[0] = trace_P

    print("iteration {:2d}: trace(P) = {:e}, trace(P)-Ne/2 = {: e}".format(
      i+1, trace_P, trace_P-options.Ne/2.0))

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
  else:
    print("failed to converge")
    sys.exit(1)

##############################################

if __name__ == "__main__":
  main()