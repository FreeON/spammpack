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
#include "timer.h"

#include <assert.h>

/** The constructor. The multiply is initialized so that a call to
 * Multiply::multiply() will result in:
 *
 * @f[ C \leftarrow C + A \times B @f]
 *
 * @param A Matrix A.
 * @param B Matrix B.
 * @param C Matrix C.
 */
Multiply::Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C)
{
  this->A = A;
  this->B = B;
  this->C = C;

  depth = -1;
  convolution = NULL;
#ifdef PRUNE_CONVOLUTION
  convolutionMap = NULL;
#endif
  PEMap = NULL;
  PEMap_norm = NULL;
  complexity = -1;
}

/** The destructor.
 */
Multiply::~Multiply (void)
{
  DEBUG("destructor\n");

  if(PEMap != NULL)
  {
    delete[] PEMap;
  }

  if(PEMap_norm != NULL)
  {
    delete[] PEMap_norm;
  }

#ifdef PRUNE_CONVOLUTION
  if(convolutionMap != NULL)
  {
    for(int tier = 0; tier <= depth; tier++)
    {
      delete[] convolutionMap[tier];
    }
    delete[] convolutionMap;
  }
#endif
}

/** Initialize the convolution space.
 *
 * @param initialPE The initial PE of the MultiplyElement chares.
 * @param alignPEs Align PEs in the diagonal matrix case.
 * @param cb The callback to send to when done.
 */
void Multiply::init (int initialPE, bool alignPEs, CkCallback &cb)
{
  MatrixNodeMsg *ANodes = A.getNodes();
  MatrixNodeMsg *BNodes = B.getNodes();
  MatrixNodeMsg *CNodes = C.getNodes();

  MatrixInfoMsg *CInfo = C.info();

  depth = CInfo->depth;

  convolution = new CProxy_MultiplyElement[depth+1];
  memset(convolution, 0, sizeof(CProxy_MultiplyElement)*(depth+1));

#ifdef PRUNE_CONVOLUTION
  convolutionMap = new bool*[depth+1];
#endif

  for(int tier = depth; tier >= 0; tier--)
  {
    int NTier = 1 << tier;

    size_t bytes = NTier*NTier*NTier
      *(sizeof(MultiplyElement)
          +CInfo->blocksize*CInfo->blocksize*sizeof(double));
    INFO("created %dx%dx%d convolution, %llu MultiplyElements "
        "using %llu bytes (%s)\n",
        NTier, NTier, NTier, (size_t) NTier*NTier*NTier,
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
          convolution[tier](i, j, k).insert(CInfo->N, CInfo->blocksize,
              CInfo->N_basic, tier, depth, ANodes->nodes[tier],
              BNodes->nodes[tier], CNodes->nodes[tier], initialPE);
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

  //for(int tier = 0; tier <= depth; tier++)
  //{
  //  convolution[tier].init(CkCallbackResumeThread());
  //}

  delete ANodes;
  delete BNodes;
  delete CNodes;

  cb.send();
}

/** Multiply two Matrix objects.
 *
 * A full convolution (as a space filling curve) is constructed, i.e. the
 * product
 *
 * @f[ C \leftarrow \beta C + \alpha A \times B @f]
 *
 * is calculated.
 *
 * @param tolerance The SpAMM tolerance.
 * @param alpha The factor @f$ \alpha @f$.
 * @param beta The factor @f$ \beta @f$.
 * @param cb The callback.
 */
void Multiply::multiply (double tolerance, double alpha, double beta, CkCallback &cb)
{
  DEBUG("tolerance = %e\n", tolerance);

  /* Prune convolution to reduce the complexity in the symbolic part of the
   * multiply. */
  MatrixNodeMsg *ANodes = A.getNodes();
  MatrixNodeMsg *BNodes = B.getNodes();

  //Timer tPrune("pruning");
  //tPrune.start();
  for(int tier = 0; tier < depth; tier++)
  {
    DEBUG("pruning tier %d\n", tier+1);
#ifdef PRUNE_CONVOLUTION
    int NTier = 1 << (tier+1);
    convolution[tier].pruneProduct(NTier, convolutionMap[tier+1], tolerance,
        ANodes->nodes[tier+1], BNodes->nodes[tier+1],
        CkCallbackResumeThread());
#else
    convolution[tier].pruneProduct(tolerance, ANodes->nodes[tier+1],
        BNodes->nodes[tier+1], CkCallbackResumeThread());
#endif
#ifdef PRUNE_CONVOLUTION
    /* Mark the next tier as complete so that the load balancer can work. */
    DEBUG("calling doneInserting() on tier %d\n", tier+1);
    convolution[tier+1].doneInserting();
#else
    DEBUG("done with tier %d\n", tier+1);
#endif
  }
  //tPrune.stop();
  //INFO("%s\n", tPrune.to_str());

  delete ANodes;
  delete BNodes;

  /* Multiply. */
  DEBUG("multiply\n");
  //Timer tMultiply("multiply");
  //tMultiply.start();
  convolution[depth].multiply(tolerance, CkCallbackResumeThread());
  //tMultiply.stop();
  //INFO("%s\n", tMultiply.to_str());

  DEBUG("scale by beta (%e)\n", beta);
  //Timer tScale("scale");
  //tScale.start();
  C.scale(beta, CkCallbackResumeThread());
  //tScale.stop();
  //INFO("%s\n", tScale.to_str());

  DEBUG("storeBack\n");
  //Timer tStoreBack("storeBack");
  //tStoreBack.start();
  convolution[depth].storeBack(alpha, CkCallbackResumeThread());
  //tStoreBack.stop();
  //INFO("%s\n", tStoreBack.to_str());

  /* Update norms. */
  DEBUG("update norms\n");
  //Timer tNorm("setNorm");
  //tNorm.start();
  C.setNorm(CkCallbackResumeThread());
  //tNorm.stop();
  //INFO("%s\n", tStoreBack.to_str());

  DEBUG("done\n");
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

/** Update the complexity count for this product.
 *
 * @param cb The callback to signal when done.
 */
void Multiply::updateComplexity (CkCallback &cb)
{
  this->cb = cb;
  CkCallback done(CkReductionTarget(Multiply, doneComplexity), thisProxy);
  convolution[depth].updateComplexity(done);
}

/** The reduction target for Multiply::updateComplexity.
 *
 * @param complexity The complexity of the Multiply.
 */
void Multiply::doneComplexity (double complexity)
{
  this->complexity = complexity;
  cb.send();
}

/** Return the complexity.
 *
 * @return The complexity.
 */
DoubleMsg * Multiply::getComplexity (void)
{
  return new DoubleMsg(complexity);
}

#include "multiply.def.h"
