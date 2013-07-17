#ifndef __MESSAGES_H
#define __MESSAGES_H

#include "messages.decl.h"

#include "node.h"

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

    CProxy_Node *root;
    MatrixInfoMsg (CProxy_Node *root);
};

#endif
