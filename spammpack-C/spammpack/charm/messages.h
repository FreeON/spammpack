#ifndef __MESSAGES_
#define __MESSAGES_

#include "messages.decl.h"
#include "matrix.h"

class EmptyMsg : public CMessage_EmptyMsg
{
  public:

    EmptyMsg ();
};

class IntMsg : public CMessage_IntMsg
{
  public:

    float i;

    IntMsg ();
    IntMsg (int i);
};

class FloatMsg : public CMessage_FloatMsg
{
  public:

    float a;

    FloatMsg ();
    FloatMsg (float a);
};

class MatrixInfoMsg : public CMessage_MatrixInfoMsg
{
  public:

    int N;
    int chunksize;
    CProxy_MatrixNode *root;

    MatrixInfoMsg (int N, int chunksize, CProxy_MatrixNode *root);
};

class NodeInfoMsg : public CMessage_NodeInfoMsg
{
  public:

    int NLower[2];
    int NUpper[2];
    int chunksize;
    float *block;
    CProxy_MatrixNode *child[4];

    NodeInfoMsg (int NLower[2], int NUpper[2], int chunksize, float *block, CProxy_MatrixNode *child[4]);
};

class BlockMsg : public CMessage_BlockMsg
{
  public:

    int chunksize;
    float *block;

    BlockMsg (int chunksize, float *block);
};

#endif
