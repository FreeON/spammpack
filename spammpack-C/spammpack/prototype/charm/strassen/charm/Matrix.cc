#include "Matrix.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <sstream>

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
    printf("index out of bounds\n");
    exit(1);
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
  if(i < 0 || j < 0 || i >= N || j >= N)
  {
    printf("index out of bounds\n");
    exit(1);
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
  MatrixMsg *A_msg = A.info();
  MatrixMsg *B_msg = B.info();

  if(A_msg->N != B_msg->N || A_msg->N != N)
  {
    printf("dimension mismatch (A = %d, B = %d, C = %d)\n", A_msg->N, B_msg->N, N);
    exit(1);
  }

  if(A_msg->blocksize != B_msg->blocksize || A_msg->blocksize != blocksize)
  {
    printf("blocksize mismatch (A = %d, B = %d, C = %d)\n",
        A_msg->blocksize, B_msg->blocksize, blocksize);
    exit(1);
  }

  if(A_msg->root == NULL || B_msg->root == NULL)
  {
    /* [FIXME] Delete C. */
    return new EmptyMsg();
  }

  if(root == NULL)
  {
    root = new CProxy_Node;
    *root = CProxy_Node::ckNew(blocksize, 0, 0, NPadded, NPadded);
  }

  return root->matmul(*A_msg->root, *B_msg->root);
}

#include "Matrix.def.h"
