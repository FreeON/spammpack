## @file
#
# Functions related to spectral bounds of matrices.
#
# @author Nicolas Bock <nicolas.bock@freeon.org>
# @author Matt Challacombe <matt.challacombe@freeon.org>

import numpy as np

## Calculate the spectral bounds using the Gershgorin circle theorem.
#
# @param A The matrix.
#
# @return The bounds.
def gershgorin (A):
  bounds = np.zeros(2, dtype = 'float')
  bounds[0] = A[0, 0]
  bounds[1] = A[0, 0]

  for i in range(A.shape[0]):
    R = 0
    for j in range(A.shape[1]):
      if i != j:
        R += abs(A[i, j])

    if A[i, i]-R < bounds[0]:
      bounds[0] = A[i, i]-R
    if A[i, i]+R > bounds[1]:
      bounds[1] = A[i, i]+R

  return bounds

## Calculate the spectral bounds using a full eigensolve.
#
# @param A The matrix.
#
# @return The bounds.
def eigensolve (A):
  (eigenvalues, _) = np.linalg.eig(A)
  bounds = np.array([ np.min(eigenvalues), np.max(eigenvalues) ])
  return bounds

## Calculate the spectral bounds.
#
# @param A The matrix.
# @param estimate_type The type to use.
#
# @return The bounds.
def spectral_bounds (A, estimate_type = gershgorin):
  if A.ndim != 2:
    raise Exception("ndim has to be 2")
  if A.shape[0] != A.shape[1]:
    raise Exception("only square matrices are supported")

  return estimate_type(A)
