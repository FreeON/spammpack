#include "node.h"

Node::Node (int depth, int childsize, int tier)
{
  this->depth = depth;
  this->childsize = childsize;
  this->tier = tier;
}

#include "node.def.h"
