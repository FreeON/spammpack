#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"

class Node : public CBase_Node
{
  private:

    int depth;
    int childsize;
    int tier;
    CkCallback normCB;
    CProxy_Node *child;

  public:

    Node (int depth, int childsize, int tier);
    Node (CkMigrateMessage *m);
    void norm (const CkCallback &cb);
    void normDone (double norm);
};

#endif
