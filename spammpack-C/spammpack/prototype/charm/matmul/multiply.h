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

    /** The tree depth of the matrix. */
    int depth;

    /** The convolution. There is a convolution for each tier. */
    CProxy_MultiplyElement *convolution;

    /** A callback. */
    CkCallback cb;

  public:

    Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
        int blocksize, int depth, CProxy_Node ANodes, CProxy_Node BNodes,
        CProxy_Node CNodes);
    ~Multiply (void);
    void multiply (double tolerance, CkCallback &cb);
    void updatePEMap (CkCallback &cb);
    void donePEMap (CkReductionMsg *data);
};

#endif
