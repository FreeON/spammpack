#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"
#include "types.h"
#include <map>
#include <list>

class Node : public CBase_Node
{
  private:

    int depth;
    int blocksize;
    int tier;

    int iLower, iUpper;
    int jLower, jUpper;

    std::map<int, bool> callbackSet;
    std::map<int, CkCallback> cb;

    std::map<int, std::list<int> > childWorking;
    CProxy_Node *child[4];

    double *block;

  public:

    Node (int depth, int blocksize, int tier,
        int iLower, int iUpper, int jLower, int jUpper);
    DoubleMsg * get (int i, int j);
    void initialize (int initType, int index, CkCallback &cb);
    void initializeDone (IntMsg *index);
    void matmul (int index, CProxy_Node A, CProxy_Node B, CkCallback &cb);
};

#endif
