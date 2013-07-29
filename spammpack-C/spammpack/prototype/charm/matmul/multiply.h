/** @file
 *
 * The header file for the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MULTIPLY_H
#define __MULTIPLY_H

#include "multiplyelement.h"

#include "multiply.decl.h"

/** A multiplication. */
class Multiply : public CBase_Multiply
{
  private:

    /** The convolution. A 3D space filling curve in the product space. */
    CProxy_MultiplyElement convolution;

    /** The callback. */
    CkCallback cb;

  public:

    Multiply ();
    void multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
        CkCallback &cb);
    void multiplyDone ();
    void storeBackDone ();
};

#endif
