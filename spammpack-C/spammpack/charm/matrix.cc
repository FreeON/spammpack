#include "matrix.h"
#include "messages.h"
#include "utilities.h"

/** Initialize matrix with random elements.
 *
 * @param N The size of the matrix.
 * @param chunksize The size of the submatrix chunks.
 */
Matrix::Matrix (int N,
    int chunksize)
{
  this->N = N;
  this->chunksize = chunksize;

  root = NULL;

  int depth = (int) floor(log(N/(double) chunksize)/log(2.));
  while(chunksize*pow(2, depth) < N) { depth++; }
  NPadded = chunksize*pow(2, depth);

  LOG_INFO("Initializing %dx%d matrix, with %dx%d chunks, "
      "depth = %d, NPadded = %d\n", N, N, chunksize, chunksize, depth,
      NPadded);
}

void Matrix::get (int i,
    int j,
    CkFuture f)
{
  GetMsg *aij = new GetMsg();

  if(this->root == NULL)
  {
    aij->a = 0;
  }

  else
  {
    CkFuture f_child = CkCreateFuture();
    root->get(i, j, f_child);
    GetMsg *result = (GetMsg*) CkWaitFuture(f_child);
    aij->a = result->a;
    delete result;
    CkReleaseFuture(f_child);
  }

  CkSendToFuture(f, aij);
}

void Matrix::set (int i,
    int j,
    float aij,
    CkFuture f)
{
  LOG_INFO("setting A(%d,%d) <- %f\n", i, j, aij);

  if(root == NULL)
  {
    /* Create root node. */
    root = new CProxy_MatrixNode;
    int NLower[2] = { 0, 0 };
    int NUpper[2] = { NPadded, NPadded };
    *root = CProxy_MatrixNode::ckNew(NLower, NUpper, chunksize);
  }

  CkFuture f_child = CkCreateFuture();
  LOG_INFO("created future\n");
  root->set(i, j, aij, f_child);
  EmptyMsg *m = (EmptyMsg*) CkWaitFuture(f_child);
  CkReleaseFuture(f_child);

  EmptyMsg *result = new EmptyMsg();
  CkSendToFuture(f_child, result);
}

/** Print a matrix. */
void Matrix::print (CkFuture f)
{
  CkFuture f_child = CkCreateFuture();
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      float aij;
      get(i, j, f_child);
      GetMsg *result = (GetMsg*) CkWaitFuture(f_child);
      CkPrintf(" %1.2f", aij);
    }
    CkPrintf("\n");
  }
}

#include "matrix.def.h"
