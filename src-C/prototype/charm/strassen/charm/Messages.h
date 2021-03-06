#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "Messages.decl.h"
#include "Node.h"

class DataMsg : public CMessage_DataMsg
{
  public:

    int blocksize;
    double *data;

    DataMsg (int blocksize, double *data);
};

class DoubleMsg : public CMessage_DoubleMsg
{
  public :

    double x;

    DoubleMsg ();
    DoubleMsg (double x);
};

class EmptyMsg : public CMessage_EmptyMsg
{
};

class IntMsg : public CMessage_IntMsg
{
  public:

    int i;

    IntMsg ();
    IntMsg (int i);
};

class MatrixMsg : public CMessage_MatrixMsg
{
  public:

    int N;
    int NPadded;
    int blocksize;
    int depth;
    CProxy_Node *root;

    MatrixMsg (int N, int NPadded, int blocksize, int depth, CProxy_Node *root);
};

class NodeMsg : public CMessage_NodeMsg
{
  public:

    int iLower;
    int iUpper;
    int jLower;
    int jUpper;
    int blocksize;
    CProxy_Node *child[4];
    double *data;

    NodeMsg (int iLower, int iUpper, int jLower, int jUpper, int blocksize,
        CProxy_Node *child[4], double *data);
};

#endif
