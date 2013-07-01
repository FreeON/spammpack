#ifndef __NODE_H
#define __NODE_H

#include "Node.decl.h"
#include "Messages.h"

class Node : public CBase_Node
{
  private:

    int blocksize;

    int iLower;
    int iUpper;
    int jLower;
    int jUpper;

    CProxy_Node *child[4];
    double *data;

  public:

    Node (int blocksize, int iLower, int jLower, int iUpper, int jUpper);
    EmptyMsg * set (int i, int j, double aij);
    DoubleMsg * get (int i, int j);
    EmptyMsg * matmul (CProxy_Node A, CProxy_Node B);
    NodeMsg * info ();
    DataMsg * getData ();
};

#endif
