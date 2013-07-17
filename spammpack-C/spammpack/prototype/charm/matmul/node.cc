#include "node.h"
#include "logger.h"
#include "messages.h"

Node::Node (int depth, int blocksize, int tier,
    int iLower, int iUpper,
    int jLower, int jUpper)
{
  this->depth = depth;
  this->blocksize = blocksize;
  this->tier = tier;
  for(int i = 0; i < 4; i++)
  {
    child[i] = NULL;
  }
  callbackSet = false;
  this->block = NULL;
}

void Node::random (int index, CkCallback &cb)
{
  LOG("(%d) random\n", index);

  if(callbackSet)
  {
    LOG("(%d) callback already set\n", index);
    CkExit();
  }

  callbackSet = true;
  this->cb = cb;

  if(tier == depth)
  {
    LOG("(%d) reached depth\n", index);
    if(block == NULL)
    {
      block = new double[blocksize*blocksize];
    }
    for(int i = 0; i < blocksize*blocksize; i++)
    {
      block[i] = rand();
    }
    cb.send(new IntMsg(index));
  }

  else
  {
    for(int i = 0; i < 4; i++)
    {
      childDone[i] = false;
      child[i] = new CProxy_Node;
      *child[i] = CProxy_Node::ckNew(depth, blocksize, tier+1,
          iLower+(iUpper-iLower)/2*(i >> 1), iLower+(iUpper-iLower)/2*((i >> 1)+1),
          jLower+(jUpper-jLower)/2*(i & 1), jLower+(jUpper-jLower)/2*((i & 1)+1));
      CkCallback thisCB(CkIndex_Node::randomDone(NULL), thisProxy);
      child[i]->random((index << 2) | i, thisCB);
    }
  }
}

void Node::randomDone (IntMsg *index)
{
  LOG("(%d) received index %d\n", index->i >> 2, index->i);
  childDone[index->i & 3] = true;
  for(int i = 0; i < 4; i++)
  {
    if(!childDone[i])
    {
      return;
    }
  }
  LOG("(%d) sending to callback\n", index->i >> 2);
  callbackSet = false;
  cb.send(new IntMsg(index->i >> 2));
}

void Node::zero (int index, CkCallback &cb)
{
  LOG("(%d) zero\n", index);

  if(callbackSet)
  {
    LOG("(%d) callback already set\n", index);
    CkExit();
  }

  callbackSet = true;
  this->cb = cb;

  if(tier == depth)
  {
    LOG("(%d) reached depth\n", index);
    if(block == NULL)
    {
      block = new double[blocksize*blocksize];
    }
    memset(block, 0, blocksize*blocksize*sizeof(double));
    cb.send(new IntMsg(index));
  }

  else
  {
    for(int i = 0; i < 4; i++)
    {
      childDone[i] = false;
      child[i] = new CProxy_Node;
      *child[i] = CProxy_Node::ckNew(depth, blocksize, tier+1,
          iLower+(iUpper-iLower)/2*(i >> 1), iLower+(iUpper-iLower)/2*((i >> 1)+1),
          jLower+(jUpper-jLower)/2*(i & 1), jLower+(jUpper-jLower)/2*((i & 1)+1));
      CkCallback thisCB(CkIndex_Node::zeroDone(NULL), thisProxy);
      child[i]->zero((index << 2) | i, thisCB);
    }
  }
}

void Node::zeroDone (IntMsg *index)
{
  LOG("(%d) received index %d\n", index->i >> 2, index->i);
  childDone[index->i & 3] = true;
  for(int i = 0; i < 4; i++)
  {
    if(!childDone[i])
    {
      return;
    }
  }
  LOG("(%d) sending to callback\n", index->i >> 2);
  callbackSet = false;
  cb.send(new IntMsg(index->i >> 2));
}

DoubleMsg * Node::get (int i, int j)
{
  if(tier == depth)
  {
  }

  else
  {
  }

  return new DoubleMsg(0);
}

#include "node.def.h"
