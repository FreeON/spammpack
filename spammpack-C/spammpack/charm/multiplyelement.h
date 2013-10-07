/** @file
 *
 * The header file for the MultiplyElement class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __MULTIPLYELEMENT_H
#define __MULTIPLYELEMENT_H

#include "config.h"
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

    /** The norm product of this element. */
    double norm_product;

    /** The result matrix. */
    double *CResult;

#ifndef PRUNE_CONVOLUTION
    /** A flag indicating whether this MultiplyElement is enabled or not (as a
     * hack while real pruning doesn't work. */
    bool isEnabled;
#endif

    /** The next tier convolution. If this MultiplyElement is on tier ==
     * depth, then the proxy is not defined. */
    CProxy_MultiplyElement nextConvolution;

  public:

    MultiplyElement (int blocksize, int tier, int depth, CProxy_Node A,
        CProxy_Node B, CProxy_Node C);
    MultiplyElement (CkMigrateMessage *msg);
    ~MultiplyElement ();
    void pup (PUP::er &p);
    void multiply (double tolerance, CkCallback &cb);
#ifdef PRUNE_CONVOLUTION
    void pruneProduct (int NTier,
        bool *nextConvolutionMap,
        double tolerance,
        CProxy_Node ANodes,
        CProxy_Node BNodes,
        CkCallback &cb);
#else
    void pruneProduct (double tolerance,
        CProxy_Node ANodes,
        CProxy_Node BNodes,
        CkCallback &cb);
#endif
    void setNextConvolution (CProxy_MultiplyElement nextConvolution);
#ifndef PRUNE_CONVOLUTION
    void enable (CkCallback &cb);
    void disable (CkCallback &cb);
#endif
    void storeBack (double alpha, CkCallback &cb);
    void PEMap (CkCallback &cb);
    void complexity (CkCallback &cb);
};

#endif
