/* @file
 *
 * The implementation of the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "matrix.h"
#include "messages.h"
#include "logger.h"

/** The constructor. */
Matrix::Matrix (int N, int blocksize)
{
  this->N = N;
  this->blocksize = blocksize;

  /* Calculate tree depth. */
  depth = -1;
  for(int i = N/blocksize; i > 0; i >>= 1)
  {
    depth++;
  }
  if(blocksize*(1 << depth) < N) depth++;
  NPadded = blocksize*(1 << depth);

  int NTier = 1 << depth;

  DEBUG("N = %d, blocksize = %d, depth = %d, NPadded = %d, NTier = %d\n",
      N, blocksize, depth, NPadded, NTier);

  nodes = CProxy_Node::ckNew(N, depth, blocksize, depth, NTier, NTier);
}

/** Get some basic information on the matrix.
 *
 * @param tier The tier to return.
 *
 * @return The matrix information.
 */
MatrixInfoMsg * Matrix::info (void)
{
  return new MatrixInfoMsg (N, blocksize, depth, NPadded, nodes);
}

/** Convert a Matrix to a dense matrix.
 *
 * @return The dense matrix.
 */
DenseMatrixMsg * Matrix::toDense (void)
{
  DenseMatrixMsg *A = new (N*N) DenseMatrixMsg();

  int NTier = 1 << depth;

  for(int i = 0; i < NTier; i++) {
    for(int j = 0; j < NTier; j++)
    {
      DenseMatrixMsg *block = nodes(i, j).getBlock();

      for(int l = i*blocksize; l < (i+1)*blocksize && l < N; l++) {
        for(int m = j*blocksize; m < (j+1)*blocksize && m < N; m++)
        {
          A->A[BLOCK_INDEX(l, m, 0, 0, N)] = block->A[BLOCK_INDEX(l, m, i*blocksize, j*blocksize, blocksize)];
        }
      }

      delete block;
    }
  }

  return A;
}

#include "matrix.def.h"
