#include "Messages.h"

DataMsg::DataMsg (int blocksize, double *data)
{
  this->blocksize = blocksize;
  memcpy(this->data, data, blocksize*blocksize*sizeof(double));
}

DoubleMsg::DoubleMsg ()
{
  x = 0;
}

DoubleMsg::DoubleMsg (double x)
{
  this->x = x;
}

IntMsg::IntMsg ()
{
  i = 0;
}

IntMsg::IntMsg (int i)
{
  this->i = i;
}

MatrixMsg::MatrixMsg (int N, int NPadded, int blocksize, int depth,
    CProxy_Node *root)
{
  this->N = N;
  this->NPadded = NPadded;
  this->blocksize = blocksize;
  this->depth = depth;
  this->root = root;
}

NodeMsg::NodeMsg (int iLower, int iUpper, int jLower, int jUpper,
    int blocksize, CProxy_Node *child[4], double *data)
{
  this->iLower = iLower;
  this->iUpper = iUpper;
  this->jLower = jLower;
  this->jUpper = jUpper;
  this->blocksize = blocksize;
  for(int i = 0; i < 4; i++)
  {
    this->child[i] = child[i];
  }
  this->data = data;
}

#include "Messages.def.h"
