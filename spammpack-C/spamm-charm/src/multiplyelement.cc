/** @file
 *
 * The implementation of the MultiplyElement class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "multiplyelement.h"

#include "chunk.h"
#include "index.h"
#include "lapack_interface.h"
#include "logger.h"
#include "messages.h"
#include "types.h"
#include "utilities.h"

#include <assert.h>
#include <bitset>

/** Some convenience macros for logging. Wrap the logging format string with
 * LB and LE. */
#define LB "tier %d ME(%d,%d,%d) "

/** Some convenience macros for logging. Wrap the logging format string with
 * LB and LE. */
#define LE , tier, thisIndex.x, thisIndex.y, thisIndex.z

/** The constructor.
 *
 * @param The matrix size.
 * @param blocksize The blocksize.
 * @param tier The tier.
 * @param depth The depth of the matrix.
 * @param A The Nodes of this tier in A.
 * @param B The Nodes of this tier in B.
 * @param C The Nodes of this tier in C.
 */
MultiplyElement::MultiplyElement (int N, int blocksize, int N_basic,
    int tier, int depth, CProxy_Node A, CProxy_Node B, CProxy_Node C)
{
  DEBUG(LB"constructor\n"LE);

  this->blocksize = blocksize;
  this->N_basic = N_basic;
  this->tier = tier;
  this->depth = depth;
  this->A = A;
  this->B = B;
  this->C = C;
  this->norm_product = 0;

  this->N = N;

  iLower = thisIndex.x*blocksize;
  jLower = thisIndex.y*blocksize;

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

  chunksize = 0;
  CResult = NULL;

#ifndef PRUNE_CONVOLUTION
  isEnabled = true;
#endif
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
MultiplyElement::MultiplyElement (CkMigrateMessage *msg)
{
  DEBUG("ME(%d,%d,%d) migration constructor\n", thisIndex.x, thisIndex.y,
      thisIndex.z);
  CResult = NULL;
}

/** The destructor.
 */
MultiplyElement::~MultiplyElement ()
{
  DEBUG(LB"destructor\n"LE);
  if(CResult != NULL)
  {
    DEBUG("free'ing chunk at %p\n", CResult);
    free(CResult);
    CResult = NULL;
    chunksize = 0;
  }
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
  p|N_basic;
  p|depth;
  p|tier;
  p|A;
  p|B;
  p|C;
  p|norm_product;
#ifndef PRUNE_CONVOLUTION
  p|isEnabled;
#endif
  p|nextConvolution;
  p|chunksize;
  p|N;
  p|iLower;
  p|jLower;
  bool resultSet = (CResult != NULL);
  p|resultSet;
  if(resultSet)
  {
    if(p.isUnpacking())
    {
      CResult = malloc(chunk_sizeof(blocksize, N_basic));
      assert(chunksize == chunk_sizeof(blocksize, N_basic));
    }
    PUParray(p, (char*) CResult, chunksize);
  }
  DEBUG(LB"pup()\n"LE);
}

/** Initialize this MultiplyElement.
 *
 * @param cb The callback to reduce to.
 */
void MultiplyElement::init (CkCallback &cb)
{
  INFO(LB"initializing\n"LE);
  contribute(cb);
}

/** Multiply nodes.
 *
 * @f[ C \leftarrow \beta C + \alpha A \times B @f]
 *
 * @param tolerance The multiplication tolerance.
 * @param cb The callback.
 */
void MultiplyElement::multiply (double tolerance, CkCallback &cb)
{
  assert(tier == depth);

  DEBUG(LB"multiply\n"LE);

#ifndef PRUNE_CONVOLUTION
  if(isEnabled)
#endif
  {
    NodeInfoMsg *AInfo = A(thisIndex.x, thisIndex.z).info();
    NodeInfoMsg *BInfo = B(thisIndex.z, thisIndex.y).info();

    norm_product = AInfo->norm*BInfo->norm;

    DEBUG(LB"multiplying blocks, ANorm*BNorm = %e\n"LE, norm_product);

    delete AInfo;
    delete BInfo;

    if(norm_product > tolerance)
    {
      if(CResult != NULL)
      {
        ABORT(LB"result is not NULL\n"LE);
      }

      DEBUG(LB"creating result chunk, N = %d\n"LE, N);
      CResult = chunk_alloc(blocksize, N_basic, N, iLower, jLower);
      chunksize = chunk_sizeof(blocksize, N_basic);

      /* Calculate C_{ij} = A_{ik} B_{kj}. */
      DEBUG(LB"requesting ChunkMsg from A and B\n"LE);
      ChunkMsg *AChunk = A(thisIndex.x, thisIndex.z).getChunk();
      ChunkMsg *BChunk = B(thisIndex.z, thisIndex.y).getChunk();

#ifdef DEBUG_OUTPUT
      chunk_print(AChunk->chunk, "tier %d ME(%d,%d,%d) AChunk(%d,%d):", tier,
          thisIndex.x, thisIndex.y, thisIndex.z, thisIndex.x, thisIndex.z);
      chunk_print(BChunk->chunk, "tier %d ME(%d,%d,%d) BChunk(%d,%d):", tier,
          thisIndex.x, thisIndex.y, thisIndex.z, thisIndex.z, thisIndex.y);
#endif

      DEBUG(LB"calling multiply on result\n"LE);
      chunk_multiply(AChunk->chunk, BChunk->chunk, CResult);

#ifdef DEBUG_OUTPUT
      /** For debugging. */
      chunk_print(CResult, LB"result:"LE);
#endif

      DEBUG(LB"deleting ChunkMsg from A and B\n"LE);
      delete AChunk;
      delete BChunk;
    }

#ifndef PRUNE_CONVOLUTION
    else
    {
      INFO(LB"this product was enabled but its norm product was below the tolerance\n"LE);
    }
#endif
  }

#ifndef PRUNE_CONVOLUTION
  else
  {
    DEBUG(LB"skipping disabled element\n"LE);
    norm_product = 0;
  }
#endif

  contribute(cb);
}

#ifdef PRUNE_CONVOLUTION
/** Prune the next tier based on the Node norms by applying the SpAMM
 * tolerance.
 *
 * @param NTier The size of the convolution array.
 * @param nextConvolutionMap A bool map indicating which MultiplyElement is
 * active.
 * @param tolerance The SpAMM tolerance.
 * @param ANodes The @link Node nodes @endlink of Matrix A on the tier below this one.
 * @param BNodes The @link Node nodes @endlink of Matrix B on the tier below this one.
 * @param cb The callback to reduce to.
 */
void MultiplyElement::pruneProduct (int NTier,
    bool *nextConvolutionMap,
    double tolerance,
    CProxy_Node ANodes,
    CProxy_Node BNodes,
    CkCallback &cb)
#else
/** Prune the next tier based on the Node norms by applying the SpAMM
 * tolerance.
 *
 * @param tolerance The SpAMM tolerance.
 * @param ANodes The @link Node nodes @endlink of Matrix A on the tier below this one.
 * @param BNodes The @link Node nodes @endlink of Matrix B on the tier below this one.
 * @param cb The callback to reduce to.
 */
void MultiplyElement::pruneProduct (double tolerance,
    CProxy_Node ANodes,
    CProxy_Node BNodes,
    CkCallback &cb)
#endif
{
  DEBUG(LB"pruning next tier\n"LE);

#ifndef PRUNE_CONVOLUTION
  if(isEnabled)
#endif
  {
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
            DEBUG(LB"keeping tier %d, convolution(%d,%d,%d) "
                "A[%d:%d,%d:%d]*B[%d:%d,%d:%d] "
                "(%e * %e = %e > %e)\n"LE, tier+1,
                nextX, nextY, nextZ,
                AInfo->iLower, AInfo->iUpper, AInfo->jLower, AInfo->jUpper,
                BInfo->iLower, BInfo->iUpper, BInfo->jLower, BInfo->jUpper,
                AInfo->norm, BInfo->norm, AInfo->norm*BInfo->norm, tolerance);
#ifdef PRUNE_CONVOLUTION
            if(!nextConvolutionMap[BLOCK_INDEX_3(nextX, nextY, nextZ, NTier)])
            {
              nextConvolution(nextX, nextY, nextZ).insert(blocksize, tier+1, depth, A, B, C);
              nextConvolutionMap[BLOCK_INDEX_3(nextX, nextY, nextZ, NTier)] = false;
            }
#else
            nextConvolution(nextX, nextY, nextZ).enable(CkCallbackResumeThread());
#endif
          }

          else
          {
            /* If necessary, destroy MultiplyElement. */
            DEBUG(LB"pruning tier %d, convolution(%d,%d,%d) "
                "A[%d:%d,%d:%d]*B[%d:%d,%d:%d] "
                "(%e * %e = %e <= %e)\n"LE, tier+1,
                nextX, nextY, nextZ,
                AInfo->iLower, AInfo->iUpper, AInfo->jLower, AInfo->jUpper,
                BInfo->iLower, BInfo->iUpper, BInfo->jLower, BInfo->jUpper,
                AInfo->norm, BInfo->norm, AInfo->norm*BInfo->norm, tolerance);
#ifdef PRUNE_CONVOLUTION
            if(nextConvolutionMap[BLOCK_INDEX_3(nextX, nextY, nextZ, NTier)])
            {
              nextConvolution(nextX, nextY, nextZ).ckDestroy();
              nextConvolutionMap[BLOCK_INDEX_3(nextX, nextY, nextZ, NTier)] = true;
            }
#else
            nextConvolution(nextX, nextY, nextZ).disable(CkCallbackResumeThread());
#endif
          }

          delete AInfo;
          delete BInfo;
        }
      }
    }
  }

  contribute(cb);
}

/** Set the convolution of the next tier.
 *
 * @param nextConvolution The chare array of @link MultiplyElement
 * MultiplyEleemnts @endlink of the next tier.
 */
void MultiplyElement::setNextConvolution (CProxy_MultiplyElement nextConvolution)
{
  this->nextConvolution = nextConvolution;
}

#ifndef PRUNE_CONVOLUTION
/** Enable this MultiplyElement.
 *
 * @param cb The callback to send to.
 */
void MultiplyElement::enable (CkCallback &cb)
{
  if(!isEnabled)
  {
    isEnabled = true;
    DEBUG(LB"enabling previously disabled element\n"LE);
  }
  cb.send();
}

/** Disable this MultiplyElement.
 *
 * @param cb The callback to send to.
 */
void MultiplyElement::disable (CkCallback &cb)
{
  if(isEnabled == true)
  {
    DEBUG(LB"disabling\n"LE);

    if(tier < depth)
    {
      /* Disable the elements below this one. */
      for(int i = 0; i < 2; i++)
      {
        int nextX = (thisIndex.x << 1) | i;
        for(int j = 0; j < 2; j++)
        {
          int nextY = (thisIndex.y << 1) | j;
          for(int k = 0; k < 2; k++)
          {
            int nextZ = (thisIndex.z << 1) | k;

            DEBUG(LB"disabling nextConvolution(%d,%d,%d)\n"LE, nextX, nextY,
                nextZ);
            nextConvolution(nextX, nextY, nextZ).disable(CkCallbackResumeThread());
          }
        }
      }
    }

    /* Disable this MultiplyElement. */
    isEnabled = false;
  }
  cb.send();
}
#endif

/** Push the C submatrices back into the C Matrix.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param cb The callback.
 */
void MultiplyElement::storeBack (double alpha, CkCallback &cb)
{
  assert(tier == depth);

#ifndef PRUNE_CONVOLUTION
  if(isEnabled)
#endif
  {
    DEBUG(LB"storing in C\n"LE);

    if(CResult != NULL)
    {
      DEBUG(LB"calling blockAdd with Chunk at %p\n"LE, CResult);
      C(thisIndex.x, thisIndex.y).chunkAdd(alpha, chunksize, (char*) CResult);

      /* Reset result for possible next iteration. */
      DEBUG(LB"deleting CResult at %p for next iteration\n"LE, CResult);
      free(CResult);
      CResult = NULL;
    }
  }

  contribute(cb);
}

/** Create a PE map of the Matrix @link Node nodes @endlink.
 *
 * @param cb The callback for the reduction.
 */
void MultiplyElement::PEMap (CkCallback &cb)
{
  DEBUG(LB"PE %d, norm = %e\n"LE, CkMyPe(), norm_product);

  struct PEMap_MultiplyElement_t *result = (struct PEMap_MultiplyElement_t*) malloc(sizeof(PEMap_MultiplyElement_t));

  result->index[0] = thisIndex.x;
  result->index[1] = thisIndex.y;
  result->index[2] = thisIndex.z;
  result->PE = CkMyPe();
  result->norm_product = norm_product;

  contribute(sizeof(PEMap_MultiplyElement_t), result, CkReduction::set, cb);
}

/** Update the complexity of this Multiply.
 *
 * @param cb The callback for the reduction.
 */
void MultiplyElement::complexity (CkCallback &cb)
{
  int complexity = 1;
#ifndef PRUNE_CONVOLUTION
  if(!isEnabled) { complexity = 0; }
#endif
  contribute(sizeof(int), &complexity, CkReduction::sum_int, cb);
}

#include "multiplyelement.def.h"
