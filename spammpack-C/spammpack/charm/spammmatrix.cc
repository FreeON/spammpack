#include "spammmatrix.h"

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

  root->set(i, j, aij);
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

MatrixNode::MatrixNode (const int NLower[2], const int NUpper[2],
    const int chunksize)
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
void MatrixNode::get (const int i, const int j, float &aij)
{
  int N = NUpper[0]-NLower[0];

  if(N == chunksize)
  {
    if(block == NULL)
    {
      aij = 0;
    }

    else
    {
      int index = (i-NLower[0])+(j-NLower[1])*chunksize;
      aij = block[index];
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
      aij = 0;
    }

    else
    {
      /* Recurse. */
      child[childindex]->get(i, j, aij);
    }
  }
}

/** Set a matrix element.
 *
 * @param i The row index.
 * @param j The column index.
 * @param aij The matrix element.
 */
void MatrixNode::set (const int i, const int j, const float aij)
{
  int N = this->NUpper[0]-this->NLower[0];

  CkPrintf("setting A[%d][%d] to %f\n", i, j, aij);

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
    child[childindex]->set(i, j, aij);
  }
}

#include "spammmatrix.def.h"
