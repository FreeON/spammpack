/** @file
 *
 * The header file for the Multiply class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MULTIPLY_H
#define __MULTIPLY_H

#include "multiply.decl.h"

/** A multiplication. */
class Multiply : public CBase_Multiply
{
  private:

    /** Matrix A. */
    CProxy_Matrix A;

    /** Matrix B. */
    CProxy_Matrix B;

    /** Matrix C. */
    CProxy_Matrix C;

    /** The convolution. */
    CProxy_MultiplyElement convolution;

  public:

    Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
        int blocksize, int depth, CProxy_Node ANodes, CProxy_Node BNodes,
        CProxy_Node CNodes);
    void multiply (double tolerance, CkCallback &cb);
};

#endif
