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

class MultiplyElement : public CBase_MultiplyElement
{
  private:

    /** The convolution index, a linear 3D index. */
    int index;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    CProxyElement_Node ANode;
    CProxyElement_Node BNode;
    CProxyElement_Node CNode;

  public:

    MultiplyElement (int blocksize,
        CProxyElement_Node ANode,
        CProxyElement_Node BNode,
        CProxyElement_Node CNode);
    MultiplyElement (CkMigrateMessage *msg);
    void multiply (CkCallback &cb);
};

class Multiply : public CBase_Multiply
{
  private:

    CProxy_MultiplyElement convolution;
    CkCallback done;

  public:

    Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C);
    void multiply (CkCallback &cb);
    void multiplyDone ();
};

#endif
