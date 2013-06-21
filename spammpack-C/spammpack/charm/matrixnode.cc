#include "matrixnode.h"
#include "messages.h"
#include "utilities.h"

/** Constructor. */
MatrixNode::MatrixNode (int NLower[2],
    int NUpper[2],
    int chunksize)
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
void MatrixNode::get (int i,
    int j,
    CkFuture f)
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

    else if (i >= NLower[0]+N/2 && j < NLower[1]+N/2)
    {
      childindex = 1;
    }

    else if(i < NLower[0]+N/2 && j >= NLower[1]+N/2)
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
void MatrixNode::set (int i,
    int j,
    float aij,
    CkFuture f)
{
  int N = this->NUpper[0]-this->NLower[0];

  LOG_INFO("setting A[%d][%d] to %f\n", i, j, aij);

  if(N == chunksize)
  {
    if(block == NULL)
    {
      block = (float*) calloc(chunksize*chunksize, sizeof(float));
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

    else if (i >= this->NLower[0]+N/2 && j < this->NLower[1]+N/2)
    {
      childindex = 1;
      NLower[0] = this->NLower[0]+N/2;
      NUpper[0] = this->NUpper[0];
      NLower[1] = this->NLower[1];
      NUpper[1] = this->NLower[1]+N/2;
    }

    else if(i < this->NLower[0]+N/2 && j >= this->NLower[1]+N/2)
    {
      childindex = 2;
      NLower[0] = this->NLower[0];
      NUpper[0] = this->NLower[0]+N/2;
      NLower[1] = this->NLower[1]+N/2;
      NUpper[1] = this->NUpper[1];
    }

    else if(i >= this->NLower[0]+N/2 && j >= this->NLower[1]+N/2)
    {
      childindex = 3;
      NLower[0] = this->NLower[0]+N/2;
      NUpper[0] = this->NUpper[0];
      NLower[1] = this->NLower[1]+N/2;
      NUpper[1] = this->NUpper[1];
    }

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

#include "matrixnode.def.h"
