#include "algebra.h"
#include "matrixnode.h"
#include "messages.h"
#include "utilities.h"

/** Constructor. */
MatrixNode::MatrixNode (int NLower[2], int NUpper[2], int chunksize)
{
  for(int i = 0; i < 2; i++)
  {
    this->NLower[i] = NLower[i];
    this->NUpper[i] = NUpper[i];
  }
  this->chunksize = chunksize;

  block = NULL;
  for(int i = 0; i < 4; i++)
  {
    child[i] = NULL;
  }
}

/** Get a matrix element.
 */
void MatrixNode::get (int i, int j, CkFuture f)
{
  int N = NUpper[0]-NLower[0];
  FloatMsg *aij = new FloatMsg();

  if(N == chunksize)
  {
    if(block == NULL)
    {
      aij->a = 0;
    }

    else
    {
      int index = (i-NLower[0])+(j-NLower[1])*chunksize;
      aij->a = block[index];
    }
  }

  else
  {
    int childindex = -1;

    if(i < NLower[0]+N/2 && j < NLower[1]+N/2)
    {
      childindex = 0;
    }

    else if(i < NLower[0]+N/2 && j >= NLower[1]+N/2)
    {
      childindex = 1;
    }

    else if (i >= NLower[0]+N/2 && j < NLower[1]+N/2)
    {
      childindex = 2;
    }

    else if(i >= NLower[0]+N/2 && j >= NLower[1]+N/2)
    {
      childindex = 3;
    }

    if(child[childindex] == NULL)
    {
      aij->a = 0;
    }

    else
    {
      /* Recurse. */
      CkFuture f_child = CkCreateFuture();
      child[childindex]->get(i, j, f_child);
      FloatMsg *result = (FloatMsg*) CkWaitFuture(f_child);
      aij->a = result->a;
      delete result;
      CkReleaseFuture(f_child);
    }
  }

  CkSendToFuture(f, aij);
}

/** Set a matrix element.
 *
 * @param i The row index.
 * @param j The column index.
 * @param aij The matrix element.
 */
void MatrixNode::set (int i, int j, float aij, CkFuture f)
{
  int N = this->NUpper[0]-this->NLower[0];

  LOG_DEBUG("N = %d (chunksize = %d)\n", N, chunksize);
  LOG_DEBUG("setting A[%d][%d] to %f\n", i, j, aij);

  if(N == chunksize)
  {
    if(block == NULL)
    {
      LOG_DEBUG("creating new block\n");
      block = new float[chunksize*chunksize];

      for(int i = 0; i < chunksize*chunksize; i++)
      {
        block[i] = 0;
      }
    }

    /* Figure out index in submatrix block. */
    int index = (i-NLower[0])+(j-NLower[1])*chunksize;
    block[index] = aij;
  }

  else
  {
    int childindex = -1;
    int NLower[2];
    int NUpper[2];

    if(i < this->NLower[0]+N/2 && j < this->NLower[1]+N/2)
    {
      childindex = 0;
      NLower[0] = this->NLower[0];
      NUpper[0] = this->NLower[0]+N/2;
      NLower[1] = this->NLower[1];
      NUpper[1] = this->NLower[1]+N/2;
    }

    else if(i < this->NLower[0]+N/2 && j >= this->NLower[1]+N/2)
    {
      childindex = 1;
      NLower[0] = this->NLower[0];
      NUpper[0] = this->NLower[0]+N/2;
      NLower[1] = this->NLower[1]+N/2;
      NUpper[1] = this->NUpper[1];
    }

    else if (i >= this->NLower[0]+N/2 && j < this->NLower[1]+N/2)
    {
      childindex = 2;
      NLower[0] = this->NLower[0]+N/2;
      NUpper[0] = this->NUpper[0];
      NLower[1] = this->NLower[1];
      NUpper[1] = this->NLower[1]+N/2;
    }

    else if(i >= this->NLower[0]+N/2 && j >= this->NLower[1]+N/2)
    {
      childindex = 3;
      NLower[0] = this->NLower[0]+N/2;
      NUpper[0] = this->NUpper[0];
      NLower[1] = this->NLower[1]+N/2;
      NUpper[1] = this->NUpper[1];
    }

    LOG_DEBUG("recursing to NLower = { %d, %d }, NUpper = { %d, %d }\n",
        NLower[0], NLower[1], NUpper[0], NUpper[1]);

    if(child[childindex] == NULL)
    {
      child[childindex] = new CProxy_MatrixNode;
      *child[childindex] = CProxy_MatrixNode::ckNew(NLower, NUpper, chunksize);
    }

    /* Recurse. */
    CkFuture f_child = CkCreateFuture();
    child[childindex]->set(i, j, aij, f_child);
    EmptyMsg *m = (EmptyMsg*) CkWaitFuture(f_child);
    CkReleaseFuture(f_child);
    delete m;
  }

  EmptyMsg *m = new EmptyMsg();
  CkSendToFuture(f, m);
}

/** Get some basic information on this node. */
NodeInfoMsg * MatrixNode::getInfo ()
{
  NodeInfoMsg *m = new NodeInfoMsg(NLower, NUpper, chunksize, block, child);

  return m;
}

/** Get the matrix blocks. */
BlockMsg * MatrixNode::getBlock ()
{
  BlockMsg *m = new (chunksize*chunksize) BlockMsg(chunksize, block);
  if(block != NULL)
  {
    for(int i = 0; i < chunksize*chunksize; i++)
    {
      LOG_DEBUG("copying block[%d] = %f\n", i, block[i]);
      m->block[i] = block[i];
    }
  }

  return m;
}

void MatrixNode::multiply (CProxy_MatrixNode A, CProxy_MatrixNode B, CkFuture f)
{
  /* Get info on nodes. */
  NodeInfoMsg *A_info = A.getInfo();
  NodeInfoMsg *B_info = B.getInfo();

  LOG_DEBUG("A: NLower = { %d, %d }, NUpper = { %d, %d }\n",
      A_info->NLower[0], A_info->NLower[1],
      A_info->NUpper[0], A_info->NUpper[1]);

  LOG_DEBUG("C: NLower = { %d, %d }, NUpper = { %d, %d }\n",
      NLower[0], NLower[1],
      NUpper[0], NUpper[1]);

  int N = NUpper[0]-NLower[0];

  LOG_DEBUG("N = %d (chunksize = %d)\n", N, chunksize);

  /* Multiply. */
  if(N == chunksize)
  {
    LOG_DEBUG("N == chunksize\n");
    if(A_info->block != NULL && B_info->block != NULL)
    {
      if(block == NULL)
      {
        LOG_DEBUG("creating new block in C\n");
        block = new float[chunksize*chunksize];
        for(int i = 0; i < chunksize*chunksize; i++)
        {
          block[i] = 0;
        }
      }
      LOG_DEBUG("receiving matrix blocks for A and B\n");
      BlockMsg *A_block = A.getBlock();
      BlockMsg *B_block = B.getBlock();

      logging::printBlock(logging::DEBUG, chunksize, "A_block", A_block->block);
      logging::printBlock(logging::DEBUG, chunksize, "B_block", B_block->block);
      logging::printBlock(logging::DEBUG, chunksize, "C_block", block);

      LOG_DEBUG("calling matmul()\n");
      matmul(chunksize, A_block->block, B_block->block, block);
    }
  }

  else
  {
    CkFuture child_f[8];
    int childIndex = 0;

    LOG_DEBUG("recursing\n");

    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++) {
        for(int k = 0; k < 2; k++)
        {
          if(A_info->child[i*2+k] != NULL && B_info->child[k*2+j] != NULL)
          {
            if(child[i*2+j] == NULL)
            {
              child[i*2+j] = new CProxy_MatrixNode();
              int NLower[2] = { this->NLower[0]+i*N/2, this->NLower[1]+j*N/2 };
              int NUpper[2] = { NLower[0]+N/2, NLower[1]+N/2 };
              *child[i*2+j] = CProxy_MatrixNode::ckNew(NLower, NUpper, chunksize);

              LOG_DEBUG("creating new C node, NLower = { %d, %d }, NUpper = { %d, %d }\n",
                  NLower[0], NLower[1], NUpper[0], NUpper[1]);
            }

            LOG_DEBUG("Creating future[%d]\n", childIndex);
            child_f[childIndex] = CkCreateFuture();
            child[i*2+j]->multiply(*A_info->child[i*2+k], *B_info->child[k*2+j], child_f[childIndex]);
            childIndex++;
          }
        }
      }
    }

    LOG_DEBUG("waiting on %d futures\n", childIndex);
    for(int i = 0; i < childIndex; i++)
    {
      EmptyMsg *m = (EmptyMsg*) CkWaitFuture(child_f[i]); delete m;
      CkReleaseFuture(child_f[i]);
      LOG_DEBUG("child_f[%d] finished\n", i);
    }
  }

  /* Return signal that we are done. */
  CkSendToFuture(f, new EmptyMsg());
}

#include "matrixnode.def.h"
