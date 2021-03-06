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

    /** The name of the matrix this node is part of. This is currently used to
     * identify the matrix. */
    char *name;

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The size of the basic submatrices. */
    int N_basic;

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

    /** The size of the chunk (necessary for message marshalling). */
    size_t chunksize;

    /** The local matrix chunk (if we are at the lowest tier). */
    void *chunk;

  public:

    Node (int N, int depth, int blocksize, int N_basic, int tier,
        int nameLength, char *name);
    Node (CkMigrateMessage *msg);
    ~Node (void);
    void pup (PUP::er &p);
    void init (CkCallback &cb);
    NodeInfoMsg * info (void);
    DenseMatrixMsg * toDense (void);
    ChunkMsg * getChunk (void);
    void set (int blocksize, int N_basic, double *A, CkCallback &cb);
    void setNorm (CProxy_Node nodes, CkCallback &cb);
    void chunkAdd (double alpha, size_t chunksize, char *chunk);
    void add (double alpha, double beta, CProxy_Node B, CkCallback &cb);
    void trace (CkCallback &cb);
    void PEMap (CkCallback &cb);
    void scale (double alpha, CkCallback &cb);
    void addIdentity (double alpha, CkCallback &cb);
};

#endif
