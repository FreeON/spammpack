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

#include "block.h"

/** The Node class. */
class Node : public CBase_Node
{
  private:

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

    /** The local matrix tree. */
    Block *block;

  public:

    Node (int N, int depth, int blocksize, int tier);
    Node (CkMigrateMessage *msg);
    ~Node (void);
    void pup (PUP::er &p);
    void init (CkCallback &cb);
    NodeInfoMsg * info (void);
    DenseMatrixMsg * toDense (void);
    void set (int blocksize, double *A, CkCallback &cb);
    void setNorm (CProxy_Node nodes, CkCallback &cb);
    void blockAdd (double alpha, int blocksize, double *A);
    void add (double alpha, double beta, CProxy_Node B, CkCallback &cb);
    void trace (CkCallback &cb);
    void PEMap (CkCallback &cb);
    void scale (double alpha, CkCallback &cb);
    void addIdentity (double alpha, CkCallback &cb);
};

#endif
