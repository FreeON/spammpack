#include "matrix.h"
#include "logger.h"
#include "types.h"

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

MatrixInfoMsg * Matrix::getInfo ()
{
  return new MatrixInfoMsg(root);
}

void Matrix::random (CkCallback &cb)
{
  LOG("generating random matrix\n");
  initialize(initRandom, cb);
}

void Matrix::zero (CkCallback &cb)
{
  LOG("setting matrix to zero\n");
  initialize(initZero, cb);
}

void Matrix::initialize (enum init_t initType, CkCallback &cb)
{
  if(root != NULL)
  {
    ERROR("root is not NULL\n");
    CkExit();
  }

  root = new CProxy_Node;
  *root = CProxy_Node::ckNew(depth, blocksize, 0, 0, NPadded, 0, NPadded);
  root->initialize(initType, 1, CkCallbackResumeThread());
  cb.send();
}

void Matrix::print (CkCallback &cb)
{
  if(root == NULL)
  {
    cb.send();
  }

  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      DoubleMsg *m = root->get(i, j);
      printf(" % 1.3f", m->x);
      delete m;
    }
    printf("\n");
  }
  cb.send();
}

void Matrix::multiply (CProxy_Matrix A, CProxy_Matrix B, CkCallback &cb)
{
  MatrixInfoMsg *AInfo = A.getInfo();
  MatrixInfoMsg *BInfo = B.getInfo();

  root->matmul(1, *AInfo->root, *BInfo->root, CkCallbackResumeThread());
  cb.send();
}

#include "matrix.def.h"
