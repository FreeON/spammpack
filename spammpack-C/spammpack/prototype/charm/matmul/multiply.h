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

    /** Submatrix A. */
    CProxyElement_Node ANode;

    /** Submatrix B. */
    CProxyElement_Node BNode;

    /** Submatrix C. */
    CProxyElement_Node CNode;

  public:

    MultiplyElement (int blocksize,
        CProxyElement_Node ANode,
        CProxyElement_Node BNode,
        CProxyElement_Node CNode);
    MultiplyElement (CkMigrateMessage *msg);
    void multiply (CkCallback &done);
};

/** A multiplication. */
class Multiply : public CBase_Multiply
{
  private:

    /** The convolution. A 3D space filling curve in the product space. */
    CProxy_MultiplyElement convolution;

    /** A callback to signal the completion of the multiplication. */
    CkCallback done;

  public:

    Multiply ();
    void multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
        CkCallback &cb);
    void multiplyDone ();
};

#endif
