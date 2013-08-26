/** @file
 *
 * The header file for the MultiplyElement class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MULTIPLYELEMENT_H
#define __MULTIPLYELEMENT_H

#include "multiplyelement.decl.h"

/** An element in the convolution curve. */
class MultiplyElement : public CBase_MultiplyElement
{
  private:

    /** The convolution index, a linear 3D index. */
    int index;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The tier of this Node. */
    int tier;

    /** Matrix A. */
    CProxy_Node A;

    /** Matrix B. */
    CProxy_Node B;

    /** Matrix C. */
    CProxy_Node C;

    /** The result matrix. */
    double *CResult;

  public:

    MultiplyElement (int blocksize, int tier, int depth, CProxy_Node A,
        CProxy_Node B, CProxy_Node C);
    MultiplyElement (CkMigrateMessage *msg);
    ~MultiplyElement ();
    void pup (PUP::er &p);
    void multiply (double tolerance, CkCallback &cb);
    void storeBack (CkCallback &cb);
    void printPE (CkCallback &cb);
};

#endif
