/** @file
 *
 * The implementation of the MultiplyElement class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "multiplyelement.h"
#include "blas_interface.h"
#include "messages.h"
#include "logger.h"
#include "index.h"

#include <bitset>

/** The constructor.
 *
 * @param blocksize The blocksize.
 * @param tier The tier.
 * @param depth The depth of the matrix.
 * @param A The node of matrix A.
 * @param B The node of matrix B.
 * @param C The node of matrix C.
 */
MultiplyElement::MultiplyElement (int blocksize, int tier, int depth,
    CProxy_Node A, CProxy_Node B, CProxy_Node C)
{
  DEBUG("tier %d ME(%d,%d,%d) constructor\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  this->blocksize = blocksize;
  this->tier = tier;
  this->depth = depth;
  this->A = A;
  this->B = B;
  this->C = C;

  /* Calculate the linear index. */
  std::bitset<8*sizeof(int)> iIndex(thisIndex.x);
  std::bitset<8*sizeof(int)> jIndex(thisIndex.y);
  std::bitset<8*sizeof(int)> kIndex(thisIndex.z);
  std::bitset<8*sizeof(int)> tempIndex(1);
  for(int i = 0; i < tier; i++)
  {
    tempIndex <<= 3;
    if(iIndex[i]) { tempIndex |= 4; }
    if(jIndex[i]) { tempIndex |= 2; }
    if(kIndex[i]) { tempIndex |= 1; }
  }
  index = tempIndex.to_ulong();

  CResult = NULL;
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  DEBUG("ME(%d,%d,%d) migration constructor\n", thisIndex.x, thisIndex.y,
      thisIndex.z);
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  DEBUG("tier %d ME(%d,%d,%d) destructor\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  delete[] CResult;
}

/** The PUP method.
 *
 * @param p The PUP::er object.
 */
void MultiplyElement::pup (PUP::er &p)
{
  CBase_MultiplyElement::pup(p);

  p|index;
  p|blocksize;
  p|depth;
  p|tier;
  p|A;
  p|B;
  p|C;

  int numberElements = (CResult == NULL ? 0 : blocksize*blocksize);
  p|numberElements;

  if(numberElements > 0)
  {
    if(p.isUnpacking())
    {
      CResult = new double[numberElements];
    }
    PUParray(p, CResult, numberElements);
  }
  else
  {
    if(p.isUnpacking()) { CResult = NULL; }
  }

  DEBUG("tier %d ME(%d,%d,%d) pup()\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);
}

/** Multiply nodes.
 *
 * @param tolerance The multiplication tolerance.
 * @param cb The callback.
 */
void MultiplyElement::multiply (double tolerance, CkCallback &cb)
{
  DEBUG("tier %d ME(%d,%d,%d) multiply\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  NodeInfoMsg *AInfo = A(thisIndex.x, thisIndex.z).info();
  NodeInfoMsg *BInfo = B(thisIndex.z, thisIndex.y).info();

  DEBUG("tier %d ME(%d,%d,%d) multiplying blocks\n", tier, thisIndex.x,
      thisIndex.y, thisIndex.z);

  if(AInfo->norm*BInfo->norm > tolerance)
  {
    if(CResult != NULL)
    {
      ABORT("tier %d ME(%d,%d,%d) CResult is not NULL\n", tier, thisIndex.x,
          thisIndex.y, thisIndex.z);
    }

    CResult = new double[blocksize*blocksize];
    memset(CResult, 0, sizeof(double)*blocksize*blocksize);

    /* Calculate C_{ij} = A_{ik} B_{kj}. */
    DenseMatrixMsg *ABlock = A(thisIndex.x, thisIndex.z).getBlock();
    DenseMatrixMsg *BBlock = B(thisIndex.z, thisIndex.y).getBlock();

#ifdef DEBUG_OUTPUT
    printDense(blocksize, ABlock->A, "tier %d ME(%d,%d,%d) ABlock(%d,%d):", tier,
        thisIndex.x, thisIndex.y, thisIndex.z, thisIndex.x, thisIndex.z);
    printDense(blocksize, BBlock->A, "tier %d ME(%d,%d,%d) BBlock(%d,%d):", tier,
        thisIndex.x, thisIndex.y, thisIndex.z, thisIndex.z, thisIndex.y);
#endif

#ifdef DGEMM
    double alpha = 1;
    double beta = 1;
    DGEMM("N", "N", &blocksize, &blocksize, &blocksize, &alpha, ABlock->A,
        &blocksize, BBlock->A, &blocksize, &beta, CResult, &blocksize);
#else
    for(int i = 0; i < blocksize; i++) {
      for(int j = 0; j < blocksize; j++) {
        for(int k = 0; k < blocksize; k++)
        {
          CResult[BLOCK_INDEX(i, j, 0, 0, blocksize)] +=
            ABlock->A[BLOCK_INDEX(i, k, 0, 0, blocksize)]
            *BBlock->A[BLOCK_INDEX(k, j, 0, 0, blocksize)];
        }
      }
    }
#endif

#ifdef DEBUG_OUTPUT
    /** For debugging. */
    printDense(blocksize, CResult, "tier %d ME(%d,%d,%d) result:", tier,
        thisIndex.x, thisIndex.y, thisIndex.z);
#endif

    delete ABlock;
    delete BBlock;
  }

  delete AInfo;
  delete BInfo;

  contribute(cb);
}

/** Prune the next tier based on the Node norms by applying the SpAMM
 * tolerance.
 *
 * @param tolerance The SpAMM tolerance.
 * @param ANodes The @link Node nodes @endlink of Matrix A on the tier below this one.
 * @param BNodes The @link Node nodes @endlink of Matrix B on the tier below this one.
 * @param NElement The number of elements in the convolution array.
 * @param convolutionExists A boolean array of size NElements^3 indication
 * whether a particular MultiplyElement exists in the convolution or not.
 * @param convolution The MultiplyElement chare array of the next tier.
 * @param cb The callback to reduce to.
 */
void MultiplyElement::pruneProduct (double tolerance,
    CProxy_Node ANodes,
    CProxy_Node BNodes,
    int NElement,
    bool *convolutionExists,
    CProxy_MultiplyElement convolution,
    CkCallback &cb)
{
  DEBUG("tier %d ME(%d,%d,%d) pruning next tier\n", tier, thisIndex.x,
      thisIndex.y, thisIndex.z);

  for(int i = 0; i < 2; i++)
  {
    int nextX = (thisIndex.x << 1) | i;
    for(int j = 0; j < 2; j++)
    {
      int nextY = (thisIndex.y << 1) | j;
      for(int k = 0; k < 2; k++)
      {
        int nextZ = (thisIndex.z << 1) | k;

        NodeInfoMsg *AInfo = ANodes(nextX, nextZ).info();
        NodeInfoMsg *BInfo = BNodes(nextZ, nextY).info();

        if(AInfo->norm*BInfo->norm > tolerance)
        {
          /* If necessary, create MultiplyElement. */
          DEBUG("tier %d ME(%d,%d,%d) keeping/creating tier %d, convolution(%d,%d,%d)\n",
              tier, thisIndex.x, thisIndex.y, thisIndex.z, tier+1,
              nextX, nextY, nextZ);
          if(!convolutionExists[nextX+nextY*NElement+nextZ*NElement*NElement])
          {
            convolution(nextX, nextY, nextZ).insert(blocksize, tier+1, depth, A, B, C);
            convolutionExists[nextX+nextY*NElement+nextZ*NElement*NElement] = true;
          }
        }

        else
        {
          /* If necessary, destroy MultiplyElement. */
          INFO("tier %d ME(%d,%d,%d) pruning tier %d, convolution(%d,%d,%d)\n",
              tier, thisIndex.x, thisIndex.y, thisIndex.z, tier+1,
              nextX, nextY, nextZ);
          if(convolutionExists[nextX+nextY*NElement+nextZ*NElement*NElement])
          {
            convolution(nextX, nextY, nextZ).ckDestroy();
            convolutionExists[nextX+nextY*NElement+nextZ*NElement*NElement] = false;
          }
        }
        delete AInfo;
        delete BInfo;
      }
    }
  }
  contribute(cb);
}

/** Push the C submatrices back into the C Matrix.
 *
 * @param cb The callback.
 */
void MultiplyElement::storeBack (CkCallback &cb)
{
  DEBUG("tier %d ME(%d,%d,%d) storing in C\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z);

  if(CResult != NULL)
  {
    C(thisIndex.x, thisIndex.y).add(blocksize, CResult);

    /* Reset result for possible next iteration. */
    delete[] CResult;
    CResult = NULL;
  }

  contribute(cb);
}

/** Print the PE this Node is on.
 *
 * @param cb The callback for the reduction.
 */
void MultiplyElement::printPE (CkCallback &cb)
{
  INFO("tier %d ME(%d,%d,%d) PE %d\n", tier, thisIndex.x, thisIndex.y,
      thisIndex.z, CkMyPe());
  contribute(cb);
}

#include "multiplyelement.def.h"
