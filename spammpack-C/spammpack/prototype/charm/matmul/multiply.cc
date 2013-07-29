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
  INFO("Multiply constructor\n");
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
  INFO("initializing multiply\n");

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

  INFO("filling %dx%dx%d chare array\n", (1 << CInfo->depth),
      (1 << CInfo->depth), (1 << CInfo->depth));

  convolution = CProxy_MultiplyElement::ckNew();
  for(int i = 0; i < (1 << CInfo->depth); i++) {
    for(int j = 0; j < (1 << CInfo->depth); j++) {
      for(int k = 0; k < (1 << CInfo->depth); k++)
      {
        INFO("inserting convolution at C(%d,%d) <- A(%d,%d) * B(%d,%d)\n",
            i, j, i, k, k, j);

        convolution(i, j, k).insert(CInfo->blocksize,
            AInfo->tierNode,
            BInfo->tierNode,
            CInfo->tierNode);

        DEBUG("inserted element ME(%d,%d,%d)\n", i, j, k);
      }
    }
  }
  convolution.doneInserting();
  INFO("done initializing convolution\n");

  /* Store callback. */
  this->cb = cb;

  INFO("multiplying\n");
  //CkCallback done(CkReductionTarget(Multiply, multiplyDone), thisProxy);
  //convolution.multiply(done);

  convolution.multiply(CkCallbackResumeThread());

  INFO("here\n");
  cb.send();
  INFO("here\n");

  delete AInfo;
  delete BInfo;
  delete CInfo;
}

/** The reduction target for the multiply method. */
void Multiply::multiplyDone ()
{
  INFO("storing result back\n");
  CkCallback done(CkReductionTarget(Multiply, storeBackDone), thisProxy);
  convolution.storeBack(done);
}

/** The reduction target for the storeBack method. */
void Multiply::storeBackDone ()
{
  INFO("done\n");
  cb.send();
}

#include "multiply.def.h"
