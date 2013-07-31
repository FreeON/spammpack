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

    /** The tree depth of the matrix. */
    int depth;

    /** The convolution. A 3D space filling curve in the product space, one
     * per tier. */
    CProxy_MultiplyElement *convolution;

  public:

    Multiply ();
    void multiply (double tolerance, CProxy_Matrix A, CProxy_Matrix B,
        CProxy_Matrix C, CkCallback &cb);
};

#endif
