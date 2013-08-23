/** @file
 *
 * The header file for the Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"

/** The Node class. */
class Node : public CBase_Node
{
  private:

    /** The random number seed. */
    unsigned int seed;

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The tier of this Node. */
    int tier;

    /** The lower row index. */
    int iLower;

    /** The upper row index. */
    int iUpper;

    /** The lower column index. */
    int jLower;

    /** The upper column index. */
    int jUpper;

    /** The linear index of this node. */
    unsigned int index;

    /** The norm of this matrix block. */
    double norm;

    /** The square of the norm of this matrix block. */
    double norm_2;

    /** The submatrix of size blocksize x blocksize. */
    double *block;

  public:

    Node (int N, int depth, int blocksize, int tier, unsigned int seed);
    Node (CkMigrateMessage *msg);
    ~Node (void);
    void pup (PUP::er &p);
    NodeInfoMsg * info (void);
    DenseMatrixMsg * getBlock (void);
    void set (int blocksize, double *A);
    void add (int blocksize, double *A);
};

#endif
