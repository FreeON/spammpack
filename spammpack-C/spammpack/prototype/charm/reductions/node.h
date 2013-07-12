#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"

class Node : public CBase_Node
{
  private:

    int depth;
    int childsize;
    int tier;
    CProxy_Node *child;

  public:

    Node (int depth, int childsize, int tier);
    Node (CkMigrateMessage *m);
};

#endif
