#include "matrix.h"

Matrix::Matrix (int depth, int childsize)
{
  root = new CProxy_Node;
  *root = CProxy_Node::ckNew(depth, childsize, 0, 1, 1);
}

void Matrix::norm ()
{
}

#include "matrix.def.h"
