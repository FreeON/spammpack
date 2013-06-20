#include "spammmatrix.h"
#include "messages.h"

/** Initialize matrix with random elements.
 *
 * @param N The size of the matrix.
 * @param chunksize The size of the submatrix chunks.
 */
SpAMMMatrix::SpAMMMatrix (const int N, const int chunksize)
{
  CkPrintf("Initializing %dx%d matrix, with %dx%d chunks\n", N, N, chunksize, chunksize);
  this->N = N;
  this->chunksize = chunksize;

  root = NULL;

  int depth = (int) floor(log(N/(double) chunksize)/log(2.));
  while(chunksize*pow(2, depth) < N) { depth++; }
  CkPrintf("depth = %d\n", depth);

  NPadded = chunksize*pow(2, depth);
  CkPrintf("NPadded = %d\n", NPadded);
}

void SpAMMMatrix::get (const int i, const int j, float &aij)
{
  if(this->root == NULL)
  {
    aij = 0;
  }

  else
  {
    root->get(i, j, aij);
  }
}

void SpAMMMatrix::set (const int i, const int j, const float aij)
{
  if(root == NULL)
  {
    /* Create root node. */
    root = new CProxy_MatrixNode;
    int NLower[2] = { 0, 0 };
    int NUpper[2] = { NPadded, NPadded };
    *root = CProxy_MatrixNode::ckNew(NLower, NUpper, chunksize);
  }

  CkFuture f = CkCreateFuture();
  root->set(i, j, aij, f);
  SetMsg *m = (SetMsg*) CkWaitFuture(f);
}

/** Print a matrix. */
void SpAMMMatrix::print ()
{
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
      float aij;
      get(i, j, aij);
      CkPrintf(" %1.2f", aij);
    }
    CkPrintf("\n");
  }
}

#include "spammmatrix.def.h"
