/* @file
 *
 * The implementation of the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

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

  tierNode.resize(depth+1);
  for(int tier = 0; tier <= depth; tier++)
  {
    int NTier = 1 << tier;
    int width = NPadded >> tier;

    tierNode[tier] = CProxy_Node::ckNew();
    for(int i = 0; i < NTier; i++) {
      for(int j = 0; j < NTier; j++)
      {
        tierNode[tier](i, j).insert(N, depth, blocksize, tier,
            i*width, (i+1)*width, j*width, (j+1)*width);
      }
    }
  }
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
          DoubleMsg *m = tierNode[depth](i, j).get(i_block, j_block);
          A->A[BLOCK_INDEX(i_block, j_block, 0, 0, N)] = m->x;
          delete m;
        }
      }
    }
  }

  return A;
}

/** Initialize a Matrix with random numbers. */
void Matrix::random (CkCallback &cb)
{
  DEBUG("generating random matrix\n");
  initialize(initRandom, cb);
}

/** Initialize a Matrix with zeros. */
void Matrix::zero (CkCallback &cb)
{
  DEBUG("setting matrix to zero\n");
  initialize(initZero, cb);
}

/** Initialize a Matrix.
 *
 * @param initType How to initialize the Matrix.
 * @param cb The callback.
 */
void Matrix::initialize (enum init_t initType, CkCallback &cb)
{
  for(int i = 0; i < (1 << depth); i++) {
    for(int j = 0; j < (1 << depth); j++)
    {
      tierNode[depth](i, j).initialize(initType, CkCallbackResumeThread());
    }
  }
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
          DoubleMsg *m = tierNode[depth](i, j).get(i_block, j_block);
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

/** Print the PEs the leafs sit on. */
void Matrix::printLeafPes (CkCallback &cb)
{
  for(int i = 0; i < (1 << depth); i++) {
    for(int j = 0; j < (1 << depth); j++)
    {
      tierNode[depth](i, j).printLeafPes(1, CkCallbackResumeThread());
    }
  }
  cb.send();
}

#include "matrix.def.h"
