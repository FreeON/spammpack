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

/** An element in the convolution curve. */
class MultiplyElement : public CBase_MultiplyElement
{
  private:

    /** The convolution index, a linear 3D index. */
    int index;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** Matrix A. */
    CProxy_Node A;

    /** Matrix B. */
    CProxy_Node B;

    /** Matrix C. */
    CProxy_Node C;

    /** The result matrix. */
    double *CResult;

    /** A counter, counting how many times this MultiplyElement was called. */
    int numberCalls;

  public:

    MultiplyElement (int blocksize, CProxy_Node A,
        CProxy_Node B, CProxy_Node C);
    ~MultiplyElement ();
    MultiplyElement (CkMigrateMessage *msg);
    virtual void pup (PUP::er &p);
    void multiply (CkCallback &cb);
    void storeBack (CkCallback &cb);
};

/** A multiplication. */
class Multiply : public CBase_Multiply
{
  private:

    /** The convolution. A 3D space filling curve in the product space. */
    CProxy_MultiplyElement convolution;

  public:

    Multiply ();
    void multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
        CkCallback &cb);
};

#endif
