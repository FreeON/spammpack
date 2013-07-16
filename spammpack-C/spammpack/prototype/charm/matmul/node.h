#ifndef __NODE_H
#define __NODE_H

#include "node.decl.h"

class Node : public CBase_Node
{
  private:

    int depth;
    int childsize;
    int tier;

  public:

    Node (int depth, int childsize, int tier);
};

#endif
