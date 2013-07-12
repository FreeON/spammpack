#include "matrix.h"
#include "messages.h"

Matrix::Matrix (int depth, int childsize)
{
  root = new CProxy_Node;
  *root = CProxy_Node::ckNew(depth, childsize, 0, 1, 1);
}

void Matrix::norm (const CkCallback &cb)
{
  normCB = cb;
  CkCallback cb2 = CkCallback(CkReductionTarget(Matrix, normDone), thisProxy);
  root->norm(cb2);
}

void Matrix::normDone (double norm)
{
  normCB.send(new DoubleMsg(norm));
}

#include "matrix.def.h"
