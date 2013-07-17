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
    childNull[i] = NULL;
  }
  this->block = NULL;
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

void Node::initialize (int initType, int index, CkCallback &cb)
{
  LOG("(%d) initializing\n", index);

  if(callbackSet[index])
  {
    LOG("(%d) callback already set\n", index);
    CkExit();
  }

  callbackSet[index] = true;
  this->cb[index] = cb;

  if(tier == depth)
  {
    LOG("(%d) reached depth\n", index);
    if(block == NULL)
    {
      block = new double[blocksize*blocksize];
    }
    switch(initType)
    {
      case initRandom:
        for(int i = 0; i < blocksize*blocksize; i++)
        {
          block[i] = rand();
        }
        break;

      case initZero:
        memset(block, 0, blocksize*blocksize*sizeof(double));
        break;

      default:
        ERROR("unknown initType\n");
        CkExit();
        break;
    }
    cb.send(new IntMsg(index));
  }

  else
  {
    std::list<int> tempDone;
    childWorking[index] = tempDone;
    for(int i = 0; i < 4; i++)
    {
      int childIndex = (index << 2) | i;
      childWorking[index].push_back(childIndex);
      child[i] = CProxy_Node::ckNew(depth, blocksize, tier+1,
          iLower+(iUpper-iLower)/2*(i >> 1), iLower+(iUpper-iLower)/2*((i >> 1)+1),
          jLower+(jUpper-jLower)/2*(i & 1), jLower+(jUpper-jLower)/2*((i & 1)+1));
      childNull[i] = false;
      CkCallback thisCB(CkIndex_Node::initializeDone(NULL), thisProxy);
      child[i].initialize(initType, childIndex, thisCB);
    }
  }
}

void Node::initializeDone (IntMsg *index)
{
  int thisIndex = index->i >> 2;
  LOG("(%d) received index %d\n", thisIndex, index->i);
  for(std::list<int>::iterator i = childWorking[thisIndex].begin();
      i != childWorking[thisIndex].end();
      i++)
  {
    if(*i == index->i)
    {
      childWorking[thisIndex].erase(i);
      break;
    }
  }

  if(childWorking[thisIndex].size() > 0)
  {
    return;
  }

  LOG("(%d) sending to callback\n", thisIndex);
  callbackSet[index->i] = false;
  cb[thisIndex].send(new IntMsg(thisIndex));
}

void Node::multiply (int index, CProxy_Node A, CProxy_Node B, CkCallback &cb)
{
  LOG("(%d) starting\n", index);
  //cb.send(new IntMsg(index));
}

void Node::multiplyDone (IntMsg *index)
{
}

#include "node.def.h"
