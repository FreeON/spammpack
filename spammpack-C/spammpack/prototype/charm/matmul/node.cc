#include "node.h"
#include "logger.h"
#include "messages.h"
#include "index.h"

Node::Node (int depth, int blocksize, int tier,
    int iLower, int iUpper,
    int jLower, int jUpper)
{
  this->depth = depth;
  this->blocksize = blocksize;
  this->tier = tier;
  this->iLower = iLower;
  this->iUpper = iUpper;
  this->jLower = jLower;
  this->jUpper = jUpper;
  for(int i = 0; i < 4; i++)
  {
    childNull[i] = true;
  }
  this->block = NULL;
}

NodeInfoMsg * Node::info ()
{
  NodeInfoMsg *m = new NodeInfoMsg();
  for(int i = 0; i < 4; i++)
  {
    m->childNull[i] = childNull[i];
    m->child[i] = child[i];
  }
  return m;
}

NodeBlockMsg * Node::getBlock ()
{
  NodeBlockMsg *m = new (blocksize*blocksize) NodeBlockMsg();
  memcpy(m->block, block, sizeof(double)*blocksize*blocksize);
  return m;
}

DoubleMsg * Node::get (int i, int j)
{
  DEBUG("tier %d, i = [%d,%d), j = [%d,%d), getting (%d,%d)\n",
      tier, iLower, iUpper, jLower, jUpper, i, j);
  if(tier == depth)
  {
    if(block == NULL)
    {
      return new DoubleMsg(0);
    }

    else
    {
      DEBUG("found block[%d] = %f\n", BLOCK_INDEX(i, j, iLower, jLower, blocksize),
          block[BLOCK_INDEX(i, j, iLower, jLower, blocksize)]);
      return new DoubleMsg(block[BLOCK_INDEX(i, j, iLower, jLower, blocksize)]);
    }
  }

  else
  {
    int iChild = (i < iLower+(iUpper-iLower)/2 ? 0 : 1);
    int jChild = (j < jLower+(jUpper-jLower)/2 ? 0 : 1);

    DEBUG("tier %d, calling child(%d,%d)\n", tier, iChild, jChild);
    if(childNull[CHILD_INDEX(iChild, jChild)])
    {
      DEBUG("child is NULL\n");
      return new DoubleMsg(0);
    }

    else
    {
      return child[CHILD_INDEX(iChild, jChild)].get(i, j);
    }
  }
}

void Node::initialize (int initType, int index, CkCallback &cb)
{
  DEBUG("(%d) initializing\n", index);

  if(callbackSet[index])
  {
    DEBUG("(%d) callback already set\n", index);
    CkExit();
  }

  callbackSet[index] = true;
  this->cb[index] = cb;

  if(tier == depth)
  {
    DEBUG("(%d) reached depth\n", index);
    if(block == NULL)
    {
      block = new double[blocksize*blocksize];
    }
    switch(initType)
    {
      case initRandom:
        for(int i = 0; i < blocksize*blocksize; i++)
        {
          block[i] = rand()/(double) RAND_MAX;
        }
        break;

      case initZero:
        memset(block, 0, blocksize*blocksize*sizeof(double));
        break;

      default:
        ABORT("unknown initType\n");
        break;
    }
    cb.send(new IntMsg(index));
  }

  else
  {
    std::list<int> tempDone;
    childWorking[index] = tempDone;
    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        int childIndex = (index << 2) | CHILD_INDEX(i, j);
        DEBUG("creating child[%d] with index %d\n", CHILD_INDEX(i, j), childIndex);
        childWorking[index].push_back(childIndex);
        child[CHILD_INDEX(i, j)] = CProxy_Node::ckNew(depth, blocksize,
            tier+1,
            iLower+(iUpper-iLower)/2*i, iLower+(iUpper-iLower)/2*(i+1),
            jLower+(jUpper-jLower)/2*j, jLower+(jUpper-jLower)/2*(j+1));
        childNull[CHILD_INDEX(i, j)] = false;
        CkCallback thisCB(CkIndex_Node::initializeDone(NULL), thisProxy);
        child[CHILD_INDEX(i, j)].initialize(initType, childIndex, thisCB);
      }
    }
  }
}

void Node::initializeDone (IntMsg *index)
{
  int thisIndex = index->i >> 2;
  DEBUG("(%d) received index %d\n", thisIndex, index->i);
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

  DEBUG("(%d) sending to callback\n", thisIndex);

  childWorking.erase(thisIndex);
  callbackSet[index->i] = false;
  CkCallback tempCb = cb[thisIndex];
  cb.erase(thisIndex);
  tempCb.send(new IntMsg(thisIndex));
}

void Node::multiply (int index, CProxy_Node A, CProxy_Node B, CkCallback &cb)
{
  DEBUG("(%d) starting\n", index);

  callbackSet[index] = true;
  this->cb[index] = cb;

  if(tier == depth)
  {
    DEBUG("(%d) reached depth\n", index);
    NodeBlockMsg *ABlock = A.getBlock();
    NodeBlockMsg *BBlock = B.getBlock();

    for(int i = 0; i < blocksize; i++) {
      for(int j = 0; j < blocksize; j++) {
        for(int k = 0; k < blocksize; k++)
        {
          block[BLOCK_INDEX(i, j, 0, 0, blocksize)] +=
            ABlock->block[BLOCK_INDEX(i, k, 0, 0, blocksize)]
            *BBlock->block[BLOCK_INDEX(k, j, 0, 0, blocksize)];
        }
      }
    }

    cb.send(new IntMsg(index));
  }

  else
  {
    NodeInfoMsg *AInfo = A.info();
    NodeInfoMsg *BInfo = B.info();

    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        if(childNull[CHILD_INDEX(i, j)])
        {
          ABORT("create new child\n");
        }

        for(int k = 0; k < 2; k++)
        {
          int childIndex = (index << 3) | (i << 2) | (j << 1) | k;
          childWorking[index].push_back(childIndex);
          CkCallback thisCB(CkIndex_Node::multiplyDone(NULL), thisProxy);
          DEBUG("(%d) descending C(%d,%d) <- A(%d,%d)*B(%d,%d)\n", index,
              i, j, i, k, k, j);
          child[CHILD_INDEX(i, j)].multiply(childIndex,
              AInfo->child[CHILD_INDEX(i, k)],
              BInfo->child[CHILD_INDEX(k, j)],
              thisCB);
        }
      }
    }
  }
}

void Node::multiplyDone (IntMsg *index)
{
  int thisIndex = index->i >> 3;
  DEBUG("(%d) received index %d\n", thisIndex, index->i);
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

  DEBUG("(%d) sending to callback\n", thisIndex);

  childWorking.erase(thisIndex);
  callbackSet[index->i] = false;
  CkCallback tempCb = cb[thisIndex];
  cb.erase(thisIndex);
  tempCb.send(new IntMsg(thisIndex));
}

#include "node.def.h"
