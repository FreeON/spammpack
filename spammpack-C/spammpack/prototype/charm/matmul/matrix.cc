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

  tierNode = new CProxy_Node[depth+1];
  for(int tier = 0; tier < depth+1; tier++)
  {
    tierNode[tier] = CProxy_Node::ckNew(N, depth, blocksize, tier,
        (1 << tier), (1 << tier));
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

/** Get some basic information on the matrix.
 *
 * @param tier The tier to return.
 *
 * @return The matrix information.
 */
MatrixInfoMsg * Matrix::info (int tier)
{
  DEBUG("getting matrix info (depth = %d, tier = %d)\n", depth, tier);
  MatrixInfoMsg *msg = new MatrixInfoMsg(N, blocksize, depth);
  msg->tierNode = tierNode[tier];
  return msg;
}

/** Initialize a Matrix with random numbers. */
void Matrix::random (CkCallback &cb)
{
  DEBUG("generating random matrix\n");
  initialize(initRandom, 0);
  cb.send();
}

/** Initialize a Matrix with zeros. */
void Matrix::zero (CkCallback &cb)
{
  DEBUG("setting matrix to zero\n");
  initialize(initZero, 0);
  cb.send();
}

/** Initialize a Matrix with zeros. */
void Matrix::decay (double decayConstant, CkCallback &cb)
{
  INFO("setting matrix to a matrix with decay (gamma = %f)\n", decayConstant);
  initialize(initDecay, decayConstant);
  cb.send();
}

/** Initialize a Matrix.
 *
 * @param initType How to initialize the Matrix.
 * @param decayConstant The decay constant for matrices with decay.
 */
void Matrix::initialize (int initType, double decayConstant)
{
  switch(initType)
  {
    case initZero:
    case initRandom:
      for(int i = 0; i < (1 << depth); i++) {
        for(int j = 0; j < (1 << depth); j++)
        {
          tierNode[depth](i, j).initialize(initType, CkCallbackResumeThread());
        }
      }
      break;

    case initDecay:
      {
        double *ADense = new double[N*N];
        double *ABuffer = new double[blocksize*blocksize];

        for(int i = 0; i < N; i++)
        {
          ADense[BLOCK_INDEX(i, i, 0, 0, N)] = 1+0.3*(rand()/(double) RAND_MAX - 0.5);
          for(int j = i+1; j < N; j++)
          {
            ADense[BLOCK_INDEX(i, j, 0, 0, N)] = exp(-fabs(i-j)/decayConstant)*ADense[BLOCK_INDEX(i, i, 0, 0, N)];
            ADense[BLOCK_INDEX(j, i, 0, 0, N)] = exp(-fabs(i-j)/decayConstant)*ADense[BLOCK_INDEX(i, i, 0, 0, N)];
          }
        }

#ifdef DEBUG_OUTPUT
        DEBUG("created dense matrix with decay\n");
        printDense(N, ADense);
#endif

        for(int i = 0; i < (1 << depth); i++) {
          for(int j = 0; j < (1 << depth); j++)
          {
            /* Reset buffer. */
            memset(ABuffer, 0, sizeof(double)*blocksize*blocksize);

            for(int i_block = i*blocksize; i_block < N && i_block < (i+1)*blocksize; i_block++) {
              for(int j_block = j*blocksize; j_block < N && j_block < (j+1)*blocksize; j_block++)
              {
                ABuffer[BLOCK_INDEX(i_block, j_block, i*blocksize, j*blocksize, blocksize)] =
                  ADense[BLOCK_INDEX(i_block, j_block, 0, 0, N)];
              }
            }
            tierNode[depth](i, j).set(blocksize*blocksize, ABuffer, CkCallbackResumeThread());
          }
        }

        delete[] ABuffer;
        delete[] ADense;
      }
      break;

    default:
      ABORT("unknown matrix type\n");
      break;
  }

  /* Build the upper tiers. */
  DEBUG("Building upper tiers\n");
  for(int tier = depth-1; tier >= 0; tier--)
  {
    DEBUG("tier %d, setting tierNode[%d]\n", tier, tier+1);
    tierNode[tier].setTierNode(tierNode[tier+1], CkCallbackResumeThread());
    tierNode[tier].updateNorms(CkCallbackResumeThread());
  }
}

/** Print a Matrix.
 */
void Matrix::print (CkCallback &cb)
{
  std::ostringstream o;
  o.setf(std::ios::scientific);

  if(N <= 32)
  {
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
  }

  else
  {
    INFO("matrix size too large for printing\n");
  }

  cb.send();
}

/** Print the PEs the leafs sit on.
 */
void Matrix::printLeafPes (CkCallback &cb)
{
  tierNode[depth].printLeafPes(CkCallbackResumeThread());
  cb.send();
}

#include "matrix.def.h"
