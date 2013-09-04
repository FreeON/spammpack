/** @file
 *
 * The implementation of the Matrix class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "matrix.h"
#include "messages.h"
#include "logger.h"
#include "utilities.h"
#include "index.h"

#include <assert.h>
#include <sstream>

/** The constructor.
 *
 * @param N The matrix size.
 * @param blocksize The SpAMM blocksize.
 * @param nameLength The strlen of the name.
 * @param name The matrix name.
 */
Matrix::Matrix (int N, int blocksize, int nameLength, char *name)
{
  this->name = strdup(name);
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

  nodes = new CProxy_Node[depth+1];
  for(int tier = 0; tier <= depth; tier++)
  {
    int NTier = 1 << tier;

    unsigned long bytes = NTier*NTier*(sizeof(Node)
        +blocksize*blocksize*sizeof(double));
    INFO("name = %s, N = %d, blocksize = %d, tier = %d, depth = %d, "
        "NPadded = %d, NTier = %d, creating %d Nodes using %d bytes (%s)\n",
        this->name, N, blocksize, tier, depth, NPadded,
        NTier, NTier*NTier, bytes, humanReadableSize(bytes).c_str());

    nodes[tier] = CProxy_Node::ckNew(N, depth, blocksize, tier, NTier, NTier);

    if(tier == depth)
    {
      PEMap = new int[NTier*NTier];
    }
  }
}

/** The destructor.
 */
Matrix::~Matrix (void)
{
  delete[] PEMap;
}

/** Get some basic information on the matrix.
 *
 * @return The matrix information.
 */
MatrixInfoMsg * Matrix::info (void)
{
  return new MatrixInfoMsg (N, blocksize, depth, NPadded);
}

/** Convert a Matrix to a dense matrix.
 *
 * @return The dense matrix.
 */
DenseMatrixMsg * Matrix::toDense (void)
{
  DenseMatrixMsg *A = new (N*N) DenseMatrixMsg();

  int NTier = 1 << depth;

  for(int i = 0; i < NTier; i++) {
    for(int j = 0; j < NTier; j++)
    {
      DenseMatrixMsg *block = nodes[depth](i, j).getBlock();

      for(int l = i*blocksize; l < (i+1)*blocksize && l < N; l++) {
        for(int m = j*blocksize; m < (j+1)*blocksize && m < N; m++)
        {
          A->A[BLOCK_INDEX(l, m, 0, 0, N)] = block->A[BLOCK_INDEX(l, m,
              i*blocksize, j*blocksize, blocksize)];
        }
      }

      delete block;
    }
  }

  return A;
}

/** Get the Node array on a particular tier.
 *
 * @param tier The tier.
 *
 * @return The Node array on that tier.
 */
MatrixNodeMsg * Matrix::getNodes (int tier)
{
  assert(tier >= 0 && tier <= depth);
  return new MatrixNodeMsg(nodes[tier]);
}

/** Print the PEs all @link Node Nodes @endlink are on.
 *
 * @param cb The callback to signal once all @link Node Nodes @endlink have
 * printed.
 */
void Matrix::updatePEMap (CkCallback &cb)
{
  this->cb = cb;
  CkCallback done(CkReductionTarget(Matrix, donePEMap), thisProxy);
  nodes[depth].PEMap(done);
}

/** The reduction target for Matrix::updatePEMap.
 *
 * @param msg The reduction message.
 */
void Matrix::donePEMap (CkReductionMsg *msg)
{
  int NTier = 1 << depth;
  CkReduction::setElement *current = (CkReduction::setElement*) msg->getData();
  while(current != NULL)
  {
    int *intPtr = (int*) &current->data;
    DEBUG("data = { %d, %d, %d }\n", intPtr[0], intPtr[1], intPtr[2]);
    PEMap[BLOCK_INDEX(intPtr[0], intPtr[1], 0, 0, NTier)] = intPtr[2];
    current = current->next();
  }

#ifdef ATOMIC_PE_MAP
  std::ostringstream o;
  o << "PEMap for matrix " << name << ":" << std::endl;
#else
  INFO("PEMap for matrix %s:\n", name);
#endif
  for(int i = 0; i < NTier; i++) {
    for(int j = 0; j < NTier; j++)
    {
#ifdef ATOMIC_PE_MAP
      o << "PEMap(" << i << "," << j << ") = "
        << PEMap[BLOCK_INDEX(i, j, 0, 0, NTier)] << std::endl;;
#else
      CkPrintf("PEMap(%d,%d) = %d\n", i, j, PEMap[BLOCK_INDEX(i, j, 0, 0, NTier)]);
#endif
    }
  }
#ifdef ATOMIC_PE_MAP
  o << "end of PEMap for matrix " << name << std::endl;
  INFO(o.str().c_str());
#else
  CkPrintf("end of PEMap for matrix %s\n", name);
#endif

  cb.send();
}

/** Set a matrix using a dense array.
 *
 * @param N The matrix size.
 * @param A The dense matrix.
 * @param cb The callback to signal once done.
 */
void Matrix::set (int N, double *A, CkCallback &cb)
{
  assert(this->N == N);

  DEBUG("setting matrix\n");

  /* Set the A matrix. */
  double *block = new double[blocksize*blocksize];

  for(int i = 0; i < NPadded/blocksize; i++) {
    for(int j = 0; j < NPadded/blocksize; j++)
    {
      memset(block, 0, sizeof(double)*blocksize*blocksize);

      for(int l = i*blocksize; l < (i+1)*blocksize && l < N; l++) {
        for(int m = j*blocksize; m < (j+1)*blocksize && m < N; m++)
        {
          block[BLOCK_INDEX(l, m, i*blocksize, j*blocksize, blocksize)] =
            A[BLOCK_INDEX(l, m, 0, 0, N)];
        }
      }

      nodes[depth](i, j).set(blocksize, block, CkCallbackResumeThread());
    }
  }
  delete[] block;

  /* Update norms. */
  thisProxy.setNorm(CkCallbackResumeThread());

  cb.send();
}

/** Update the norms based on the norm information of the leaf @link Node
 * nodes @endlink.
 *
 * @param cb The callback to send back to.
 */
void Matrix::setNorm (CkCallback &cb)
{
  /* Update the norms on the upper tiers. */
  for(int tier = depth-1; tier >= 0; tier--)
  {
    nodes[tier].setNorm(nodes[tier+1], CkCallbackResumeThread());
  }
  cb.send();
}

#include "matrix.def.h"
