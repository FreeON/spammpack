#include "matrix.h"
#include "messages.h"
#include "utilities.h"

/** Initialize matrix with random elements.
 *
 * @param N The size of the matrix.
 * @param chunksize The size of the submatrix chunks.
 */
Matrix::Matrix (int N, int chunksize)
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

GetMsg* Matrix::get (int i, int j)
{
  if(this->root == NULL)
  {
    return new GetMsg(0);
  }

  else
  {
    CkFuture f_child = CkCreateFuture();
    root->get(i, j, f_child);
    GetMsg *result = (GetMsg*) CkWaitFuture(f_child);
    GetMsg *m = new GetMsg(result->a);
    delete result;
    return m;
  }
}

EmptyMsg* Matrix::set (int i, int j, float aij)
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
  delete m;
  CkReleaseFuture(f_child);

  return new EmptyMsg();
}

/** Print a matrix. */
EmptyMsg* Matrix::print ()
{
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      GetMsg *m = get(i, j);
      CkPrintf(" %1.2f", m->a);
      delete m;
    }
    CkPrintf("\n");
  }
}

#include "matrix.def.h"
