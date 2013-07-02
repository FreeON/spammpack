#include "Matrix.h"
#include "Utilities.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <sstream>

Matrix::Matrix (int N, int blocksize)
{
  LOG_DEBUG("full constructor\n");

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

/** Get some information on the matrix.
 *
 * @return The MatrixMsg object.
 */
MatrixMsg * Matrix::info ()
{
  return new MatrixMsg(N, blocksize, root);
}

/** Get a matrix element.
 *
 * @param i Row index.
 * @param j Column index.
 *
 * @return The matrix element.
 */
DoubleMsg * Matrix::get (int i, int j)
{
  DoubleMsg *aij;

  if(i < 0 || j < 0 || i >= N || j >= N)
  {
    LOG_ERROR("index out of bounds\n");
    CkExit();
  }

  if(root == NULL) aij = 0;
  else aij = root->get(i, j);

  return aij;
}

/** Set a matrix element.
 *
 * @param i Row index.
 * @param j Column index.
 *
 * @return A message signifying completion.
 */
EmptyMsg * Matrix::set (int i, int j, double aij)
{
  LOG_DEBUG("setting matrix element A(%d,%d) <- %f\n", i, j, aij);

  if(i < 0 || j < 0 || i >= N || j >= N)
  {
    LOG_ERROR("index out of bounds\n");
    CkExit();
  }

  if(root == NULL)
  {
    root = new CProxy_Node;
    *root = CProxy_Node::ckNew(blocksize, 0, 0, NPadded, NPadded);
  }
  EmptyMsg *msg = root->set(i, j, aij);
  return msg;
}

/** Multiply two matrices.
 *
 * @param A Matrix A.
 * @param B Matrix B.
 *
 * @return A message signifying completion.
 */
EmptyMsg * Matrix::matmul (CProxy_Matrix A, CProxy_Matrix B)
{
  MatrixMsg *AInfo = A.info();
  MatrixMsg *BInfo = B.info();

  if(AInfo->N != BInfo->N || AInfo->N != N)
  {
    LOG_ERROR("dimension mismatch (A = %d, B = %d, C = %d)\n", AInfo->N, BInfo->N, N);
    CkExit();
  }

  if(AInfo->blocksize != BInfo->blocksize || AInfo->blocksize != blocksize)
  {
    LOG_ERROR("blocksize mismatch (A = %d, B = %d, C = %d)\n",
        AInfo->blocksize, BInfo->blocksize, blocksize);
    CkExit();
  }

  if(AInfo->root == NULL || BInfo->root == NULL)
  {
    LOG_ERROR("[FIXME] Delete C\n");
    return new EmptyMsg();
  }

  if(root == NULL)
  {
    root = new CProxy_Node;
    *root = CProxy_Node::ckNew(blocksize, 0, 0, NPadded, NPadded);
  }

  root->matmul(*AInfo->root, *BInfo->root, 0, CkCallbackResumeThread());

  return new EmptyMsg();
}

#include "Matrix.def.h"
