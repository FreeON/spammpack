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
  DEBUG("initializing multiply, tolerance = %e\n", tolerance);

  MatrixInfoMsg *AInfo = A.info();
  MatrixInfoMsg *BInfo = B.info();
  MatrixInfoMsg *CInfo = C.info();

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
    INFO("creating new convolutions\n");

    convolution = new CProxy_MultiplyElement[depth+1];

    /* First the lowest tier. */
    convolution[depth] = CProxy_MultiplyElement::ckNew(CInfo->blocksize,
        depth, depth, AInfo->tierNode[depth], BInfo->tierNode[depth],
        CInfo->tierNode[depth], 1 << depth, 1 << depth,
        1 << depth);

    /* The the upper tiers. */
    for(int i = depth-1; i >= 0; i--)
    {
      int convolutionSize = 1 << i;
      DEBUG("filling %dx%dx%d chare array\n", convolutionSize,
          convolutionSize, convolutionSize);

      convolution[i] = CProxy_MultiplyElement::ckNew(CInfo->blocksize, i,
          depth, AInfo->tierNode[i+1], BInfo->tierNode[i+1],
          CInfo->tierNode[i+1], convolutionSize, convolutionSize,
          convolutionSize);
      convolution[i].setNextTier(convolution[i+1], AInfo->tierNode[i+1],
          BInfo->tierNode[i+1], CkCallbackResumeThread());
    }
  }

  delete AInfo;
  delete BInfo;
  delete CInfo;

  DEBUG("done initializing convolution\n");

#ifdef USE_REDUCTION_TARGET
  /* Store callback. */
  this->cb = cb;

  DEBUG("multiplying\n");
  CkCallback done(CkReductionTarget(Multiply, multiplyDone), thisProxy);
  convolution.multiply(done);
#else
  for(int i = 0; i < depth+1; i++)
  {
    DEBUG("tier %d: multiplying\n", i);
    convolution[i].multiply(tolerance, CkCallbackResumeThread());
    if(i < depth)
    {
      /* In case the last reduction inserted new MultiplyElements, we need to
       * tell the load balancer. */
      DEBUG("tier %d: calling doneInserting() on tier %d\n", i, i+1);
      convolution[i+1].doneInserting();
    }
  }

  DEBUG("storing result back\n");
  convolution[depth].storeBack(CkCallbackResumeThread());

  DEBUG("sending back\n");
  cb.send();
#endif
}

#ifdef USE_REDUCTION_TARGET
/** The reduction target for the multiply method. */
void Multiply::multiplyDone ()
{
  DEBUG("storing result back\n");
  CkCallback done(CkReductionTarget(Multiply, storeBackDone), thisProxy);
  convolution[depth].storeBack(done);
}

/** The reduction target for the storeBack method. */
void Multiply::storeBackDone ()
{
  DEBUG("sending back\n");
  cb.send();
}
#endif

#include "multiply.def.h"
