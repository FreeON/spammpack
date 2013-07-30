#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "node.h"

#include "messages.decl.h"

class DenseMatrixMsg : public CMessage_DenseMatrixMsg
{
  public:

    double *A;
};

class DoubleMsg : public CMessage_DoubleMsg
{
  public:

    double x;

    DoubleMsg (double x);
};

class EmptyMsg : public CMessage_EmptyMsg
{
};

class IntMsg : public CMessage_IntMsg
{
  public:

    int i;

    IntMsg (int i);
};

class MatrixInfoMsg : public CMessage_MatrixInfoMsg
{
  public:

    /** The matrix size. */
    int N;

    /** The submatrix size at the lowest tier. */
    int blocksize;

    /** The tree depth of the matrix. */
    int depth;

    /** The array of nodes at tier == depth. */
    CProxy_Node *tierNode;

    MatrixInfoMsg (int N, int blocksize, int depth);
};

class NodeBlockMsg : public CMessage_NodeBlockMsg
{
  public:

    double *block;
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

    /** The nodes of the next tier. */
    CProxy_Node tierNode;

    NodeInfoMsg (int index, double norm, double norm_2);
};

#endif
