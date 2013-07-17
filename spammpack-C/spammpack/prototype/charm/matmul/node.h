#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"

class Node : public CBase_Node
{
  private:

    int depth;
    int blocksize;
    int tier;

    int iLower, iUpper;
    int jLower, jUpper;

    bool callbackSet;
    CkCallback cb;

    bool childDone[4];
    CProxy_Node *child[4];

    double *block;

  public:

    Node (int depth, int blocksize, int tier,
        int iLower, int iUpper, int jLower, int jUpper);
    void random (int index, CkCallback &cb);
    void randomDone (IntMsg *index);
    void zero (int index, CkCallback &cb);
    void zeroDone (IntMsg *index);
    DoubleMsg * get (int i, int j);
};

#endif
