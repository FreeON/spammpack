/** @file
 *
 * The implementation of the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "multiply.h"
#include "messages.h"
#include "logger.h"

/** The constructor.
 */
Multiply::Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
    int blocksize, int depth, CProxy_Node ANodes, CProxy_Node BNodes,
    CProxy_Node CNodes)
{
  DEBUG("Multiply constructor\n");

  this->A = A;
  this->B = B;
  this->C = C;

  int NTier = 1 << depth;

  this->convolution = CProxy_MultiplyElement::ckNew(blocksize, depth, depth,
      ANodes, BNodes, CNodes, NTier, NTier, NTier);

  DEBUG("Multiply constructor, created %dx%dx%d convolution\n", NTier, NTier,
      NTier);
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
  DEBUG("tolerance = %e\n", tolerance);
  convolution.multiply(tolerance, CkCallbackResumeThread());

  DEBUG("storeBack\n");
  convolution.storeBack(CkCallbackResumeThread());

  cb.send();
}

#include "multiply.def.h"
