#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "messages.decl.h"

#include "node.h"

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
    CProxy_Node tierNode;

    MatrixInfoMsg (int N, int blocksize, int depth, CProxy_Node tierNode);
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

    NodeInfoMsg (int index);
};

#endif
