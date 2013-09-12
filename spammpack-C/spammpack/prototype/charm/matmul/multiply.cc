/** @file
 *
 * The implementation of the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "multiply.h"
#include "multiplyelement.h"
#include "messages.h"
#include "logger.h"
#include "utilities.h"
#include "types.h"
#include "index.h"

#include <assert.h>

/** The constructor. The multiply is initialized so that a call to
 * Multiply::multiply() will result in:
 *
 * @f[ C \leftarrow C + A \times B @f]
 *
 * @param initialPE The initial PE of the MultiplyElement chares.
 * @param alignPEs Align PEs in the diagonal matrix case.
 * @param A Matrix A.
 * @param B Matrix B.
 * @param C Matrix C.
 * @param blocksize The SpAMM blocksize.
 * @param depth The depth of the matrix trees.
 * @param ANodes The Node objects of A.
 * @param BNodes The Node objects of B.
 * @param CNodes The Node objects of C.
 */
Multiply::Multiply (int initialPE, bool alignPEs, CProxy_Matrix A,
    CProxy_Matrix B, CProxy_Matrix C, int blocksize, int depth,
    CProxy_Node ANodes, CProxy_Node BNodes, CProxy_Node CNodes)
{
  this->A = A;
  this->B = B;
  this->C = C;

  this->depth = depth;

  convolution = new CProxy_MultiplyElement[depth+1];
#ifdef PRUNE_CONVOLUTION
  convolutionMap = new bool*[depth+1];
#endif
  for(int tier = depth; tier >= 0; tier--)
  {
    int NTier = 1 << tier;

    unsigned long bytes = NTier*NTier*NTier*(sizeof(MultiplyElement)+blocksize*blocksize*sizeof(double));
    INFO("created %dx%dx%d convolution, %d MultiplyElements "
        "using %d bytes (%s)\n",
        NTier, NTier, NTier, NTier*NTier*NTier,
        bytes, humanReadableSize(bytes).c_str());

    convolution[tier] = CProxy_MultiplyElement::ckNew();
#ifdef PRUNE_CONVOLUTION
    convolutionMap[tier] = new bool[NTier*NTier*NTier];
#endif
    for(int i = 0; i < NTier; i++) {
      for(int j = 0; j < NTier; j++) {
        for(int k = 0; k < NTier; k++)
        {
          if(alignPEs && i == j && i == k)
          {
            initialPE = i%CkNumPes();
          }
          convolution[tier](i, j, k).insert(blocksize, tier, depth,
              ANodes, BNodes, CNodes, initialPE);
#ifdef PRUNE_CONVOLUTION
          convolutionMap[tier][BLOCK_INDEX_3(i, j, k, NTier)] = true;
#endif
        }
      }
    }
    convolution[tier].doneInserting();

    if(tier < depth)
    {
      convolution[tier].setNextConvolution(convolution[tier+1]);
    }

    if(tier == depth)
    {
      PEMap = new int[NTier*NTier*NTier];
      PEMap_norm = new double[NTier*NTier*NTier];
    }
  }
}

/** The destructor.
 */
Multiply::~Multiply (void)
{
  delete[] PEMap;
  delete[] PEMap_norm;
#ifdef PRUNE_CONVOLUTION
  for(int tier = 0; tier <= depth; tier++)
  {
    delete[] convolutionMap[tier];
  }
  delete[] convolutionMap;
#endif
}

/** Multiply two Matrix objects.
 *
 * A full convolution (as a space filling curve) is constructed.
 *
 * @param tolerance The SpAMM tolerance.
 * @param cb The callback.
 */
void Multiply::multiply (double tolerance, CkCallback &cb)
{
  INFO("tolerance = %e\n", tolerance);

  for(int tier = 0; tier < depth; tier++)
  {
    /* Prune convolution to achieve reduced complexity in symbolic part of the
     * multiply. */
    INFO("pruning tier %d\n", tier+1);
    MatrixNodeMsg *ANodes = A.getNodes(tier+1);
    MatrixNodeMsg *BNodes = B.getNodes(tier+1);
#ifdef PRUNE_CONVOLUTION
    int NTier = 1 << (tier+1);
    convolution[tier].pruneProduct(NTier, convolutionMap[tier+1], tolerance,
        ANodes->nodes, BNodes->nodes, CkCallbackResumeThread());
#else
    convolution[tier].pruneProduct(tolerance, ANodes->nodes, BNodes->nodes,
        CkCallbackResumeThread());
#endif
    delete ANodes;
    delete BNodes;

#ifdef PRUNE_CONVOLUTION
    /* Mark the next tier as complete so that the load balancer can work. */
    INFO("calling doneInserting() on tier %d\n", tier+1);
    convolution[tier+1].doneInserting();
#endif
  }

  /* Multiply. */
  convolution[depth].multiply(tolerance, CkCallbackResumeThread());

  INFO("storeBack\n");
  convolution[depth].storeBack(CkCallbackResumeThread());

  /* Update norms. */
  INFO("update norms\n");
  C.setNorm(CkCallbackResumeThread());

  INFO("done\n");
  cb.send();
}

/** Print the PEs all @link Node Nodes @endlink are on.
 *
 * @param cb The callback to signal once all @link Node Nodes @endlink have
 * printed.
 */
void Multiply::updatePEMap (CkCallback &cb)
{
  this->cb = cb;
  CkCallback done(CkReductionTarget(Multiply, donePEMap), thisProxy);
  convolution[depth].PEMap(done);
}

/** The reduction target for Multiply::updatePEMap.
 *
 * @param msg The reduction message.
 */
void Multiply::donePEMap (CkReductionMsg *msg)
{
  int NTier = 1 << depth;

  CkReduction::setElement *current = (CkReduction::setElement*) msg->getData();
  while(current != NULL)
  {
    assert(current->dataSize == sizeof(struct PEMap_MultiplyElement_t));
    PEMap_MultiplyElement_t *result = (PEMap_MultiplyElement_t*) &current->data;
    DEBUG("data = ME(%d,%d,%d) PE %d norm %e\n", result->index[0],
        result->index[1], result->index[2], result->PE, result->norm_product);
    int matrix_offset = BLOCK_INDEX_3(result->index[0], result->index[1], result->index[2], NTier);
    PEMap[matrix_offset] = result->PE;
    PEMap_norm[matrix_offset] = result->norm_product;
    current = current->next();
  }

  cb.send();
}

/** Return the PEMap. Call udpatePEMap() first.
 *
 * @return a PEMapMsg with the PEMap.
 */
PEMapMsg * Multiply::getPEMap (void)
{
  int NTier = 1 << depth;
  PEMapMsg *msg = new (NTier*NTier*NTier, NTier*NTier*NTier) PEMapMsg();
  memcpy(msg->PEMap, PEMap, NTier*NTier*NTier*sizeof(int));
  memcpy(msg->PEMap_norm, PEMap_norm, NTier*NTier*NTier*sizeof(double));
  return msg;
}

#include "multiply.def.h"
