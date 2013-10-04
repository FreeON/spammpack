import numpy as np

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

def eigensolve (A):
  (eigenvalues, _) = np.linalg.eig(A)
  bounds = np.array([ np.min(eigenvalues), np.max(eigenvalues) ])
  return bounds

def spectral_bounds (A, estimate_type = gershgorin):
  if A.ndim != 2:
    raise Exception("ndim has to be 2")
  if A.shape[0] != A.shape[1]:
    raise Exception("only square matrices are supported")

  return estimate_type(A)
