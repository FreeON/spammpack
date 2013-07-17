#include "matrix.h"
#include "logger.h"
#include "types.h"
#include <sstream>

Matrix::Matrix (int N, int blocksize)
{
  this->N = N;
  this->blocksize = blocksize;
  this->rootNull = true;

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
  if(rootNull)
  {
    return new MatrixInfoMsg();
  }
  else
  {
    return new MatrixInfoMsg(root);
  }
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
  if(!rootNull)
  {
    ERROR("root is not NULL\n");
    CkExit();
  }

  root = CProxy_Node::ckNew(depth, blocksize, 0, 0, NPadded, 0, NPadded);
  root.initialize(initType, 1, CkCallbackResumeThread());
  cb.send();
}

void Matrix::print (CkCallback &cb)
{
  if(rootNull)
  {
    cb.send();
  }

  std::ostringstream o;
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      DoubleMsg *m = root.get(i, j);
      o << " " << m->x;
      delete m;
    }
    o << std::endl;
  }
  CkPrintf(o.str().c_str());
  cb.send();
}

void Matrix::multiply (CProxy_Matrix A, CProxy_Matrix B, CkCallback &cb)
{
  MatrixInfoMsg *AInfo = A.getInfo();
  MatrixInfoMsg *BInfo = B.getInfo();

  if(AInfo->rootNull || AInfo->rootNull)
  {
    LOG("nothing to multiply\n");
    cb.send();
  }

  LOG("starting multiply\n");
  root.multiply(1, AInfo->root, BInfo->root, CkCallbackResumeThread());
  LOG("done\n");
  cb.send();
}

#include "matrix.def.h"
