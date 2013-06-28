#include <stdlib.h>
#include "strassenOMP.h"

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

void Node::set (int i, int j, double aij)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
    }
    data[(i-iLower)*blocksize+(j-jLower)] = aij;
  }

  else
  {
    int childIndex = 0;
    int width = (iUpper-iLower)/2;
    int newILower = iLower;
    int newJLower = jLower;
    if(iLower+width <= i)
    {
      childIndex |= 1;
      newILower = iLower+width;
    }
    if(jLower+width/2 <= j)
    {
      childIndex |= 2;
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
