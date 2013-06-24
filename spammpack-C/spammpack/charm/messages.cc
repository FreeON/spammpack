#include "messages.h"

EmptyMsg::EmptyMsg () {}

IntMsg::IntMsg ()
{
  this->i = 0;
}

IntMsg::IntMsg (int i)
{
  this->i = i;
}

FloatMsg::FloatMsg ()
{
  this->a = 0;
}

FloatMsg::FloatMsg (float a)
{
  this->a = a;
}

MatrixInfoMsg::MatrixInfoMsg (int N, int chunksize, CProxy_MatrixNode *root)
{
  this->N = N;
  this->chunksize = chunksize;
  this->root = root;
}

NodeInfoMsg::NodeInfoMsg (int NLower[2], int NUpper[2], int chunksize, float *block, CProxy_MatrixNode *child[4])
{
  for(int i = 0; i < 2; i++)
  {
    this->NLower[i] = NLower[i];
    this->NUpper[i] = NUpper[i];
  }
  this->chunksize = chunksize;
  this->block = block;
  for(int i = 0; i < 4; i++)
  {
    this->child[i] = child[i];
  }
}

BlockMsg::BlockMsg (int chunksize, float *block)
{
  this->chunksize = chunksize;
  //this->block = new float[chunksize*chunksize];
}

#include "messages.def.h"
