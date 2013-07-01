#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>
#include <sstream>

#include "strassenOMP.h"

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

/** Convert a matrix in dense row-major storage to a Matrix.
 *
 * @param N The matrix dimension.
 * @param A Pointer to the dense matrix.
 */
void Matrix::convert (int N, double *A)
{
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      set(i, j, A[i*N+j]);
    }
  }
}

void Matrix::random ()
{
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      set(i, j, rand()/(double) RAND_MAX);
    }
  }
}

void Matrix::zero ()
{
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      set(i, j, 0.0);
    }
  }
}

void Matrix::set (int i, int j, double aij)
{
  if(i < 0 || j < 0 || i >= N || j >= N)
  {
    printf("index out of bounds\n");
    exit(1);
  }

  if(root == NULL)
  {
    root = new Node(blocksize, 0, 0, NPadded, NPadded);
  }
  root->set(i, j, aij);
}

double Matrix::get (int i, int j)
{
  if(i < 0 || j < 0 || i >= N || j >= N)
  {
    printf("index out of bounds\n");
    exit(1);
  }

  if(root == NULL) return 0;
  else return root->get(i, j);
}

void Matrix::print (std::string name)
{
  std::ostringstream o;
  o.setf(std::ios::fixed);
  o << name << std::endl;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      o << " " << get(i, j);
    }
    o << std::endl;
  }
  printf("%s", o.str().c_str());
}

void Matrix::matmul (Matrix A, Matrix B)
{
  if(A.N != B.N || A.N != N)
  {
    printf("dimension mismatch (A = %d, B = %d, C = %d)\n", A.N, B.N, N);
    exit(1);
  }

  if(A.blocksize != B.blocksize || A.blocksize != blocksize)
  {
    printf("blocksize mismatch (A = %d, B = %d, C = %d)\n",
        A.blocksize, B.blocksize, blocksize);
    exit(1);
  }

  if(A.root == NULL || B.root == NULL)
  {
    /* [FIXME] Delete C. */
    return;
  }

  if(root == NULL)
  {
    root = new Node(blocksize, 0, 0, NPadded, NPadded);
  }

#pragma omp parallel
  {
#pragma omp single
    {
#pragma omp task untied default(none) shared(A, B)
      root->matmul(*A.root, *B.root);
#pragma omp taskwait
    }
  }
}
