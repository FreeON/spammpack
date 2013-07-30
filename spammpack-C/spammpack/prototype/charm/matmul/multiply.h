/** @file
 *
 * The header file for the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MULTIPLY_H
#define __MULTIPLY_H

#include "config.h"

#include "multiplyelement.h"

#include "multiply.decl.h"

/** A multiplication. */
class Multiply : public CBase_Multiply
{
  private:

    /** The tree depth of the matrix. */
    int depth;

    /** The convolution. A 3D space filling curve in the product space, one
     * per tier. */
    CProxy_MultiplyElement *convolution;

#ifdef USE_REDUCTION_TARGET
    /** The callback. */
    CkCallback cb;
#endif

  public:

    Multiply ();
    void multiply (double tolerance, CProxy_Matrix A, CProxy_Matrix B,
        CProxy_Matrix C, CkCallback &cb);
#ifdef USE_REDUCTION_TARGET
    void multiplyDone ();
    void storeBackDone ();
#endif
};

#endif
