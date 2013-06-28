#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "strassenOMP.h"

Matrix::Matrix (int N, int blocksize)
{
  this->N = N;
  this->blocksize = blocksize;
  this->root = NULL;

  /* Calculate tree depth. */
  depth = -1;
  for(int i = N; i > 0; i >>= 1)
  {
    depth++;
  }
  if(1 << depth < N) depth++;
  printf("depth = %d, NPadded = %d\n", depth, 1 << depth);
  NPadded = (1 << depth);
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

  if(root == NULL)
  {
    root = new Node(blocksize, 0, 0, NPadded, NPadded);
  }
  root->set(i, j, aij);
}
