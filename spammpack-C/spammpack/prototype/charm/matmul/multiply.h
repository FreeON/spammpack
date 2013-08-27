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

    /** An array indicating whether a particular MultiplyElement exists or
     * not. Unless I find a better way, say by inquiring the proxy directly,
     * it seems that an additional array is necessary to avoid destroying an
     * already non-existant element. */
    bool **convolutionExists;

    /** The convolution. There is a convolution for each tier. */
    CProxy_MultiplyElement *convolution;

  public:

    Multiply (CProxy_Matrix A, CProxy_Matrix B, CProxy_Matrix C,
        int blocksize, int depth, CProxy_Node ANodes, CProxy_Node BNodes,
        CProxy_Node CNodes);
    ~Multiply (void);
    void multiply (double tolerance, CkCallback &cb);
    void printPE (CkCallback &cb);
};

#endif
