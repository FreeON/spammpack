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

    int tier;

    /* Store productIndex -> { childProductIndex, ... } so that we can track
     * at what point a recursion is done. */

    CkCallback parentDone;

    CProxy_Node *child[4];
    double *data;

  public:

    Node (int tier, int blocksize, int iLower, int jLower, int iUpper, int jUpper);
    NodeMsg   * info ();
    DataMsg   * getData ();
    EmptyMsg  * set (int i, int j, double aij);
    EmptyMsg  * setBlock (int iLower_arg, int jLower_arg, int iUpper_arg, int jUpper_arg, const double *ABlock);
    DoubleMsg * get (int i, int j);
#ifdef CALLBACK
    void matmul (CProxy_Node A, CProxy_Node B, int productIndex, CkCallback &done);
    void matmulDone (IntMsg *index);
#elif defined(FUTURES)
    void matmul (CProxy_Node A, CProxy_Node B, int productIndex, CkFuture f);
#else
#error "FIXME"
#endif
};

#endif
