/** @file
 *
 * The implementation of the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

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
void Multiply::multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
    CkCallback &cb)
{
  DEBUG("initializing multiply\n");

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

  if(CkMyPe() != 0)
  {
    ABORT("not on PE 0\n");
  }

  int convolutionSize = 1 << CInfo->depth;
  DEBUG("filling %dx%dx%d chare array\n", convolutionSize,
      convolutionSize, convolutionSize);

  convolution = CProxy_MultiplyElement::ckNew(CInfo->blocksize,
      AInfo->tierNode, BInfo->tierNode, CInfo->tierNode, convolutionSize,
      convolutionSize, convolutionSize);

  delete AInfo;
  delete BInfo;
  delete CInfo;

  DEBUG("done initializing convolution\n");

  /* Store callback. */
  this->cb = cb;

  DEBUG("multiplying\n");
  convolution.multiply(CkCallbackResumeThread());

  DEBUG("store back\n");
  convolution.storeBack(CkCallbackResumeThread());

  DEBUG("sending back\n");
  cb.send();
}

#include "multiply.def.h"
