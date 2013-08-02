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
#include "index.h"
#include <string.h>

/** The constructor.
 */
Multiply::Multiply ()
{
  depth = -1;
  convolution = NULL;
  DEBUG("Multiply constructor\n");
}

/** Multiply two Matrix objects.
 *
 * A full convolution (as a space filling curve) is constructed.
 *
 * @param A Matrix A.
 * @param B Matrix B.
 * @param C Matrix C.
 * @param cb The callback.
 */
void Multiply::multiply (double tolerance, CProxy_Matrix A, CProxy_Matrix B,
    CProxy_Matrix C, CkCallback &cb)
{
  DEBUG("multiply, tolerance = %e\n", tolerance);

  MatrixInfoMsg *AInfo = A.info(0);
  MatrixInfoMsg *BInfo = B.info(0);
  MatrixInfoMsg *CInfo = C.info(0);

  if(AInfo->N != BInfo->N || AInfo->N != CInfo->N)
  {
    ABORT("matrix dimension mismatch\n");
  }

  if(AInfo->blocksize != BInfo->blocksize || AInfo->blocksize != CInfo->blocksize)
  {
    ABORT("blocksize mismatch\n");
  }

  if(depth < 0)
  {
    depth = CInfo->depth;
  }

  if(depth != CInfo->depth)
  {
    ABORT("depth changed between calls to this method\n");
  }

  if(CkMyPe() != 0)
  {
    INFO("not on PE 0\n");
  }

  if(convolution == NULL)
  {
    DEBUG("creating new convolutions\n");

    convolution = new CProxy_MultiplyElement[depth+1];

    MatrixInfoMsg *AInfoTier = A.info(depth);
    MatrixInfoMsg *BInfoTier = B.info(depth);
    MatrixInfoMsg *CInfoTier = C.info(depth);

    /* First the lowest tier. */
    int convolutionSize = 1 << depth;
    DEBUG("tier %d, filling %dx%dx%d chare array\n", depth, convolutionSize,
        convolutionSize, convolutionSize);
    convolution[depth] = CProxy_MultiplyElement::ckNew(CInfo->blocksize,
        depth, depth, AInfoTier->tierNode, BInfoTier->tierNode,
        CInfoTier->tierNode, convolutionSize, convolutionSize,
        convolutionSize);

    delete AInfoTier;
    delete BInfoTier;
    delete CInfoTier;

    /* The the upper tiers. */
    for(int tier = depth-1; tier >= 0; tier--)
    {
      int convolutionSize = 1 << tier;

      DEBUG("tier %d, filling %dx%dx%d chare array\n", tier, convolutionSize,
          convolutionSize, convolutionSize);

      MatrixInfoMsg *AInfoTier = A.info(tier+1);
      MatrixInfoMsg *BInfoTier = B.info(tier+1);
      MatrixInfoMsg *CInfoTier = C.info(tier+1);

      convolution[tier] = CProxy_MultiplyElement::ckNew(CInfo->blocksize,
          tier, depth, AInfoTier->tierNode, BInfoTier->tierNode,
          CInfoTier->tierNode, convolutionSize, convolutionSize,
          convolutionSize);
      convolution[tier].setNextTier(convolution[tier+1], AInfoTier->tierNode,
          BInfoTier->tierNode, CkCallbackResumeThread());

      DEBUG("tier %d, done with convolution\n", tier);

      DEBUG("tier %d, waiting for tier to finish\n", tier);
      CkWaitQD();

      delete AInfoTier;
      delete BInfoTier;
      delete CInfoTier;
    }
  }

  delete AInfo;
  delete BInfo;
  delete CInfo;

  DEBUG("done initializing convolution\n");

  DEBUG("starting multiply\n");
  for(int tier = 0; tier < depth+1; tier++)
  {
    DEBUG("tier %d: multiplying\n", tier);
    convolution[tier].multiply(tolerance, CkCallbackResumeThread());
    if(tier < depth)
    {
      /* In case the reduction inserted or destroyed MultiplyElements, we need
       * to tell the load balancer. */
#ifdef PRUNE_CONVOLUTION
      DEBUG("tier %d: calling doneInserting() on tier %d\n", tier, tier+1);
      convolution[tier+1].doneInserting();
#endif
    }

    /* Hang out until this tier is done. */
    DEBUG("tier %d: waiting for tier to finish\n", tier);
    CkWaitQD();
  }

  DEBUG("storing result back\n");
  convolution[depth].storeBack(CkCallbackResumeThread());

  DEBUG("sending back\n");
  cb.send();
}

void Multiply::getComplexity (CkCallback &cb)
{
  DEBUG("getting complexity\n");
  this->cb = cb;
  convolution[depth].getComplexity(CkCallback(CkReductionTarget(Multiply, getComplexityDone), thisProxy));
}

void Multiply::getComplexityDone (int complexity)
{
  CkPrintf("complexity = %d block products\n", complexity);
  cb.send();
}

#include "multiply.def.h"
