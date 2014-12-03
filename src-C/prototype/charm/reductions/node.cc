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
  CkPrintf("migrating...\n");
}

void Node::norm (const CkCallback &cb)
{
  if(tier == depth)
  {
    double norm = 1.2;
    CkPrintf("returning norm %f\n", norm);
    contribute(sizeof(double), &norm, CkReduction::sum_double, cb);
  }

  else
  {
    normCB = cb;
    CkCallback cb2 = CkCallback(CkReductionTarget(Node, normDone), thisProxy);
    child->norm(cb2);
  }
}

void Node::normDone (double norm)
{
  contribute(sizeof(double), &norm, CkReduction::sum_double, normCB);
}

#include "node.def.h"
