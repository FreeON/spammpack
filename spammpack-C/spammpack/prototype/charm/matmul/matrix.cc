#include "matrix.h"

Matrix::Matrix (int N, int blocksize)
{
  this->N = N;
  this->blocksize = blocksize;
  this->root = NULL;

  /* Calculate tree depth. */
  depth = -1;
  for(int i = N/blocksize; i > 0; i >>= 1)
  {
    depth++;
  }
  if(blocksize*(1 << depth) < N) depth++;
  NPadded = blocksize*(1 << depth);
}

#include "matrix.def.h"
