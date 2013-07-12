#include "matrix.h"

Matrix::Matrix (int depth, int childsize)
{
  root = new CProxy_Node;
  *root = CProxy_Node::ckNew(depth, childsize, 0, 1, 1);
}

void Matrix::norm (const CkCallback &cb)
{
  CkPrintf("sending norm\n");
  cb.send();
}

#include "matrix.def.h"
