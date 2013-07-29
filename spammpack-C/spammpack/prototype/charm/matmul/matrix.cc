/* @file
 *
 * The implementation of the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "matrix.h"
#include "logger.h"
#include "types.h"
#include "index.h"
#include <sstream>

/** The constructor. */
Matrix::Matrix (int N, int blocksize)
{
  this->N = N;
  this->blocksize = blocksize;

  /* Calculate tree depth. */
  depth = -1;
  for(int i = N/blocksize; i > 0; i >>= 1)
  {
    depth++;
  }
  if(blocksize*(1 << depth) < N) depth++;
  NPadded = blocksize*(1 << depth);

  INFO("N = %d, blocksize = %d, depth = %d, NPadded = %d\n",
      N, blocksize, depth, NPadded);

  int NTier = 1 << depth;
  int width = NPadded >> depth;

  if(CkMyPe() != 0)
  {
    INFO("not on PE 0\n");
  }

#ifdef DENSE_ARRAYS
  tierNode = CProxy_Node::ckNew(N, depth, blocksize, depth, NTier, NTier);
#else
  tierNode = CProxy_Node::ckNew();
#endif
}

/** Convert a Matrix to a dense array. */
DenseMatrixMsg * Matrix::getDense ()
{
  DenseMatrixMsg *A = new (N*N) DenseMatrixMsg();

  memset(A->A, 0, sizeof(double)*N*N);

  for(int i = 0; i < (1 << depth); i++) {
    for(int j = 0; j < (1 << depth); j++) {
      for(int i_block = i*blocksize; i_block < (i+1)*blocksize && i_block < N; i_block++) {
        for(int j_block = j*blocksize; j_block < (j+1)*blocksize && j_block < N; j_block++)
        {
          DoubleMsg *m = tierNode(i, j).get(i_block, j_block);
          A->A[BLOCK_INDEX(i_block, j_block, 0, 0, N)] = m->x;
          delete m;
        }
      }
    }
  }

  return A;
}

/** Get some basic information on the matrix.
 *
 * @return The matrix information.
 */
MatrixInfoMsg * Matrix::info ()
{
  MatrixInfoMsg *msg = new MatrixInfoMsg(N, blocksize, depth, tierNode);
  return msg;
}

/** Initialize a Matrix with random numbers. */
void Matrix::random (CkCallback &cb)
{
  DEBUG("generating random matrix\n");
  thisProxy.initialize(initRandom, cb);
}

/** Initialize a Matrix with zeros. */
void Matrix::zero (CkCallback &cb)
{
  DEBUG("setting matrix to zero\n");
  thisProxy.initialize(initZero, cb);
}

/** Initialize a Matrix with zeros.
 *
 * @param gamma The decay constant.
 */
void Matrix::decay (double gamma, CkCallback &cb)
{
  DEBUG("setting matrix to a matrix with decay (gamma = %f)\n", gamma);
  double *ADense = new double[N*N];

  for(int i = 0; i < N; i++)
  {
    ADense[BLOCK_INDEX(i, i, 0, 0, N)] = 1+0.3*(rand()/(double) RAND_MAX - 0.5);
    for(int j = i+1; j < N; j++)
    {
      ADense[BLOCK_INDEX(i, j, 0, 0, N)] = exp(-fabs(i-j)/gamma)*ADense[BLOCK_INDEX(i, i, 0, 0, N)];
      ADense[BLOCK_INDEX(j, i, 0, 0, N)] = exp(-fabs(i-j)/gamma)*ADense[BLOCK_INDEX(i, i, 0, 0, N)];
    }
  }

  DEBUG("created dense matrix\n");
#ifdef DEBUG_OUTPUT
  for(int i = 0; i < N; i++) {
    for(int j = 0; j < N; j++)
    {
    }
  }
#endif

  double *ABuffer = new double[blocksize*blocksize];
  for(int i = 0; i < (1 << depth); i++) {
    for(int j = 0; j < (1 << depth); j++) {
      for(int i_block = i*blocksize; i_block < N && i_block < (i+1)*blocksize; i_block++) {
        for(int j_block = j*blocksize; j_block < N && j_block < (j+1)*blocksize; j_block++)
        {
          ABuffer[BLOCK_INDEX(i_block-i*blocksize, j_block-j*blocksize, 0, 0, blocksize)] =
            ADense[BLOCK_INDEX(i_block, j_block, i*blocksize, j*blocksize, N)];
#ifndef DENSE_ARRAYS
          tierNode(i, j).insert(N, depth, blocksize, depth);
#endif
          tierNode(i, j).set(blocksize*blocksize, ADense, CkCallbackResumeThread());
        }
      }
    }
  }

  delete[] ABuffer;
  delete[] ADense;
}

/** Initialize a Matrix.
 *
 * @param initType How to initialize the Matrix.
 * @param cb The callback.
 */
void Matrix::initialize (int initType, CkCallback &cb)
{
  for(int i = 0; i < (1 << depth); i++) {
    for(int j = 0; j < (1 << depth); j++)
    {
#ifndef DENSE_ARRAYS
      tierNode(i, j).insert(N, depth, blocksize, depth);
#endif
      tierNode(i, j).initialize(initType, CkCallbackResumeThread());
    }
  }
#ifndef DENSE_ARRAYS
  tierNode.doneInserting();
#endif
  cb.send();
}

/** Print a Matrix.
 */
void Matrix::print (CkCallback &cb)
{
  std::ostringstream o;
  o.setf(std::ios::fixed);

  for(int i = 0; i < (1 << depth); i++) {
    for(int i_block = i*blocksize; i_block < (i+1)*blocksize && i_block < N; i_block++) {
      for(int j = 0; j < (1 << depth); j++) {
        for(int j_block = j*blocksize; j_block < (j+1)*blocksize && j_block < N; j_block++)
        {
          DoubleMsg *m = tierNode(i, j).get(i_block, j_block);
          o << " " << m->x;
          delete m;
        }
      }
      o << std::endl;
    }
  }
  CkPrintf(o.str().c_str());
  cb.send();
}

/** Print the PEs the leafs sit on.
 */
void Matrix::printLeafPes (CkCallback &cb)
{
  tierNode.printLeafPes(CkCallbackResumeThread());
  cb.send();
}

#include "matrix.def.h"
