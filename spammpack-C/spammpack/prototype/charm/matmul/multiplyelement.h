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

    /** A counter, counting how many times this MultiplyElement was called. */
    int numberCalls;

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

    /** A mask that indicates whether the 8 convolution elements below this
     * tier exist or not. */
    bool nextConvolutionExists[2][2][2];

    /** The convolution below this tier. */
    CProxy_MultiplyElement nextConvolution;

    /** The A matrix below this tier. */
    CProxy_Node nextA;

    /** The B matrix below this tier. */
    CProxy_Node nextB;

    /** The result matrix. */
    double *CResult;

    /** A flag indicating whether this element was migrated. */
    bool wasMigrated;

  public:

    MultiplyElement (int blocksize, int tier, int depth, CProxy_Node A,
        CProxy_Node B, CProxy_Node C);
    ~MultiplyElement ();
    MultiplyElement (CkMigrateMessage *msg);
    void setNextTier (CProxy_MultiplyElement nextConvolution,
        CProxy_Node nextA, CProxy_Node nextB, CkCallback &cb);
    virtual void pup (PUP::er &p);
    void multiply (double tolerance, CkCallback &cb);
    void storeBack (CkCallback &cb);
};

#endif
