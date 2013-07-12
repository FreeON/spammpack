#include "node.h"

Node::Node (int depth, int childsize, int tier)
{
  this->depth = depth;
  this->childsize = childsize;
  this->tier = tier;

  if(tier < depth)
  {
    child = new CProxy_Node;
    *child = CProxy_Node::ckNew(depth, childsize, tier+1, childsize, childsize);
  }

  else { child = NULL; }
}

Node::Node (CkMigrateMessage *m)
{
}

#include "node.def.h"
