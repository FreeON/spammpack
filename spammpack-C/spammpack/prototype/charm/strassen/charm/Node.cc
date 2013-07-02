#include "Node.h"

#include <stdlib.h>
#include <stdio.h>
#include <cstring>

inline
int blockIndex (int i, int j, int iLower, int jLower, int blocksize)
{
  return (i-iLower)*blocksize+(j-jLower);
}

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

/** Get information on node.
 *
 * @return A NodeMsg object.
 */
NodeMsg * Node::info ()
{
  return new NodeMsg (iLower, iUpper, jLower, jUpper, blocksize, child);
}

/** Get the matrix data from a Node.
 *
 * @return A DataMsg object.
 */
DataMsg * Node::getData ()
{
  return new DataMsg();
}

DoubleMsg * Node::get (int i, int j)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL) return new DoubleMsg(0);
    else return new DoubleMsg(data[blockIndex(i, j, iLower, jLower, blocksize)]);
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
    if(child[childIndex] == NULL) return new DoubleMsg(0);
    else return child[childIndex]->get(i, j);
  }
}

EmptyMsg * Node::set (int i, int j, double aij)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
      memset(data, 0, blocksize*blocksize*sizeof(double));
    }
    data[blockIndex(i, j, iLower, jLower, blocksize)] = aij;
    return new EmptyMsg();
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
      child[childIndex] = new CProxy_Node;
      *child[childIndex] = CProxy_Node::ckNew(blocksize, newILower, newJLower,
          newILower+width, newJLower+width);
    }
    return child[childIndex]->set(i, j, aij);
  }
}

EmptyMsg * Node::matmul (CProxy_Node A, CProxy_Node B)
{
  int width = iUpper-iLower;
  NodeMsg *A_info = A.info();
  NodeMsg *B_info = B.info();

  if(width == blocksize)
  {
    DataMsg *A_data = A.getData();
    DataMsg *B_data = B.getData();

    if(A_data->data == NULL || B_data->data == NULL)
    {
      /* [FIXME] delete C. */
      return new EmptyMsg();
    }
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
      memset(data, 0, blocksize*blocksize*sizeof(double));
    }
    for(int i = iLower; i < iUpper; i++) {
      for(int j = jLower; j < jUpper; j++) {
        for(int k = A_info->jLower; k < A_info->jUpper; k++)
        {
          data[blockIndex(i, j, iLower, jLower, blocksize)] +=
            A_data->data[blockIndex(i, k, A_info->iLower, A_info->jLower, A_info->blocksize)]
            *B_data->data[blockIndex(k, j, B_info->iLower, B_info->jLower, B_info->blocksize)];
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
          if(A_info->child[childIndexA] == NULL || B_info->child[childIndexB] == NULL)
          {
            /* [FIXME] delete C. */
            continue;
          }
          if(child[childIndex] == NULL)
          {
            child[childIndex] = new CProxy_Node;
            *child[childIndex] = CProxy_Node::ckNew(blocksize,
                iLower+width/2*i, jLower+width/2*j, iLower+width/2*(i+1),
                jLower+width/2*(j+1));
          }
          child[childIndex]->matmul(*A_info->child[childIndexA], *B_info->child[childIndexB]);
        }
      }
    }
  }

  return new EmptyMsg();
}

#include "Node.def.h"
