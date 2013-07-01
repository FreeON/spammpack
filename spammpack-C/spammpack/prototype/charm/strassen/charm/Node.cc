#include "Node.h"

#include <stdlib.h>
#include <stdio.h>
#include <cstring>

Node::Node (int blocksize, int iLower, int jLower, int iUpper, int jUpper)
{
  this->blocksize = blocksize;
  this->iLower = iLower;
  this->iUpper = iUpper;
  this->jLower = jLower;
  this->jUpper = jUpper;
  for(int i = 0; i < 4; i++)
  {
    child[i] = NULL;
  }
  data = NULL;
}

inline
int Node::blockIndex (int i, int j)
{
  return (i-iLower)*blocksize+(j-jLower);
}

void Node::set (int i, int j, double aij)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
      memset(data, 0, blocksize*blocksize*sizeof(double));
    }
    data[blockIndex(i, j)] = aij;
  }

  else
  {
    int childIndex = 0;
    int width = (iUpper-iLower)/2;
    int newILower = iLower;
    int newJLower = jLower;
    if(iLower+width <= i)
    {
      childIndex |= 2;
      newILower = iLower+width;
    }
    if(jLower+width <= j)
    {
      childIndex |= 1;
      newJLower = jLower+width;
    }
    if(child[childIndex] == NULL)
    {
      child[childIndex] = new Node(blocksize, newILower, newJLower,
          newILower+width, newJLower+width);
    }
    child[childIndex]->set(i, j, aij);
  }
}

double Node::get (int i, int j)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL) return 0;
    else return data[blockIndex(i, j)];
  }

  else
  {
    int childIndex = 0;
    int width = (iUpper-iLower)/2;
    if(iLower+width <= i)
    {
      childIndex |= 2;
    }
    if(jLower+width <= j)
    {
      childIndex |= 1;
    }
    if(child[childIndex] == NULL) return 0;
    else return child[childIndex]->get(i, j);
  }
}

void Node::matmul (Node A, Node B)
{
  int width = iUpper-iLower;

  if(width == blocksize)
  {
    if(A.data == NULL || B.data == NULL)
    {
      /* [FIXME] delete C. */
      return;
    }
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
      memset(data, 0, blocksize*blocksize*sizeof(double));
    }
    for(int i = iLower; i < iUpper; i++) {
      for(int j = jLower; j < jUpper; j++) {
        for(int k = A.jLower; k < A.jUpper; k++)
        {
          data[blockIndex(i, j)] += A.data[A.blockIndex(i, k)]*B.data[B.blockIndex(k, j)];
        }
      }
    }
  }

  else
  {
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        int childIndex = (i << 1) | j;
        for(int k = 0; k < 2; k++)
        {
          int childIndexA = (i << 1) | k;
          int childIndexB = (k << 1) | j;
          if(A.child[childIndexA] == NULL || B.child[childIndexB] == NULL)
          {
            /* [FIXME] delete C. */
            continue;
          }
          if(child[childIndex] == NULL)
          {
            child[childIndex] = new Node(blocksize, iLower+width/2*i, jLower+width/2*j, iLower+width/2*(i+1), jLower+width/2*(j+1));
          }
          child[childIndex]->matmul(*A.child[childIndexA], *B.child[childIndexB]);
        }
      }
    }
  }
}

#include "Node.def.h"
