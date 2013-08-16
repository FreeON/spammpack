#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "messages.decl.h"
#include "node.h"

class MatrixInfoMsg : public CMessage_MatrixInfoMsg
{
  public:

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The padded size of the matrix. */
    int NPadded;

    /** The Node matrix at tier == depth. */
    CProxy_Node nodes;

    MatrixInfoMsg (int N, int blocksize, int depth, int NPadded, CProxy_Node nodes);
    bool equal (MatrixInfoMsg *b);
};

class NodeInfoMsg : public CMessage_NodeInfoMsg
{
  public:

    /** The linear index of this node. */
    int index;

    /** The norm of this matrix block. */
    double norm;

    /** The square of the norm of this matrix block. */
    double norm_2;

    NodeInfoMsg (int index, double norm, double norm_2);
};

class DenseMatrixMsg : public CMessage_DenseMatrixMsg
{
  public:

    double *A;
};

#endif
