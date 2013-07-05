#include "Node.h"
#include "Utilities.h"

#include <bitset>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>

inline
int blockIndex (int i, int j, int iLower, int jLower, int blocksize)
{
  return (i-iLower)*blocksize+(j-jLower);
}

Node::Node (int tier, int blocksize, int iLower, int jLower, int iUpper, int jUpper)
{
  this->blocksize = blocksize;
  this->iLower = iLower;
  this->iUpper = iUpper;
  this->jLower = jLower;
  this->jUpper = jUpper;
  for(int i = 0; i < 4; i++)
  {
    child[i] = NULL;
  }
  data = NULL;
  this->numberChildProducts = 0;
  this->tier = tier;
}

/** Get information on node.
 *
 * @return A NodeMsg object.
 */
NodeMsg * Node::info ()
{
  return new NodeMsg (iLower, iUpper, jLower, jUpper, blocksize, child);
}

/** Get the matrix data from a Node.
 *
 * @return A DataMsg object.
 */
DataMsg * Node::getData ()
{
  return new (blocksize*blocksize) DataMsg(blocksize, data);
}

DoubleMsg * Node::get (int i, int j)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL) return new DoubleMsg(0);
    else return new DoubleMsg(data[blockIndex(i, j, iLower, jLower, blocksize)]);
  }

  else
  {
    int childIndex = 0;
    int width = (iUpper-iLower)/2;
    if(iLower+width <= i)
    {
      childIndex |= 2;
    }
    if(jLower+width <= j)
    {
      childIndex |= 1;
    }
    if(child[childIndex] == NULL) return new DoubleMsg(0);
    else return child[childIndex]->get(i, j);
  }
}

EmptyMsg * Node::set (int i, int j, double aij)
{
  if(iUpper-iLower == blocksize)
  {
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
      memset(data, 0, blocksize*blocksize*sizeof(double));
    }
    data[blockIndex(i, j, iLower, jLower, blocksize)] = aij;
    return new EmptyMsg();
  }

  else
  {
    int childIndex = 0;
    int width = (iUpper-iLower)/2;
    int newILower = iLower;
    int newJLower = jLower;
    if(iLower+width <= i)
    {
      childIndex |= 2;
      newILower = iLower+width;
    }
    if(jLower+width <= j)
    {
      childIndex |= 1;
      newJLower = jLower+width;
    }
    if(child[childIndex] == NULL)
    {
      child[childIndex] = new CProxy_Node;
      *child[childIndex] = CProxy_Node::ckNew(tier+1, blocksize, newILower,
          newJLower, newILower+width, newJLower+width);
    }
    return child[childIndex]->set(i, j, aij);
  }
}

/** Multiply two matrices.
 *
 * @param A Node A.
 * @param B Node B.
 *
 * @return A message indicating completion.
 */
void Node::matmul (CProxy_Node A, CProxy_Node B, int index, CkCallback &done)
{
  int width = iUpper-iLower;
  NodeMsg *AInfo = A.info();
  NodeMsg *BInfo = B.info();

  IntMsg *indexMsg = new IntMsg(index);
  int CIndex = (index & 6) >> 1;
  int kIndex = (index & 1);
  matmulIndex = index;
  parentDone = done;

  std::bitset<20> bitIndex = std::bitset<20>(index);

  LOG_DEBUG("(%d:%d:%s) starting multiply (C index = %d), "
      "C(%d:%d,%d:%d) += A(%d:%d,%d:%d)*B(%d:%d,%d:%d)\n",
      tier, index,
      bitIndex.to_string().c_str(),
      CIndex,
      iLower+1, iUpper, jLower+1, jUpper,
      AInfo->iLower+1, AInfo->iUpper, AInfo->jLower+1, AInfo->jUpper,
      BInfo->iLower+1, BInfo->iUpper, BInfo->jLower+1, BInfo->jUpper);

  if(width == blocksize)
  {
    DataMsg *AData = A.getData();
    DataMsg *BData = B.getData();

    if(AData->data == NULL || BData->data == NULL)
    {
      LOG_ERROR("[FIXME] (%d:%d:%s) delete C, sending %d to parentDone.\n",
          tier, index, bitIndex.to_string().c_str(), index);
      parentDone.send(new EmptyMsg());
    }
    if(data == NULL)
    {
      data = new double[blocksize*blocksize];
      memset(data, 0, blocksize*blocksize*sizeof(double));
    }

    LOG_DEBUG("(%d:%d:%s) block multiply\n", tier, index, bitIndex.to_string().c_str());
    for(int i = iLower; i < iUpper; i++) {
      for(int j = jLower; j < jUpper; j++) {
        for(int k = AInfo->jLower; k < AInfo->jUpper; k++)
        {
          data[blockIndex(i, j, iLower, jLower, blocksize)] +=
            AData->data[blockIndex(i, k, AInfo->iLower, AInfo->jLower, AInfo->blocksize)]
            *BData->data[blockIndex(k, j, BInfo->iLower, BInfo->jLower, BInfo->blocksize)];
        }
      }
    }
    /* Signal that we are done. */
    LOG_DEBUG("(%d:%d:%s) sending %d to parentDone\n", tier, index, bitIndex.to_string().c_str(), index);
    parentDone.send(indexMsg);
  }

  else
  {
    LOG_DEBUG("(%d:%d:%s) descending...\n", tier, index, bitIndex.to_string().c_str());

    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        int childIndex = (i << 1) | j;
        for(int k = 0; k < 2; k++)
        {
          int childIndexA = (i << 1) | k;
          int childIndexB = (k << 1) | j;
          int productIndex = (i << 2) | (j << 1) | k;
          childProducts[numberChildProducts++] = (index << 3) | productIndex;
          if(AInfo->child[childIndexA] == NULL || BInfo->child[childIndexB] == NULL)
          {
            LOG_ERROR("(%d:%d:%s) [FIXME] delete C (product %d is NULL).\n",
                tier, index, bitIndex.to_string().c_str(), productIndex);
            thisProxy.matmulDone(new IntMsg((index << 3) | productIndex));
            continue;
          }
          if(child[childIndex] == NULL)
          {
            child[childIndex] = new CProxy_Node;
            *child[childIndex] = CProxy_Node::ckNew(tier+1, blocksize,
                iLower+width/2*i, jLower+width/2*j, iLower+width/2*(i+1),
                jLower+width/2*(j+1));
          }
          LOG_DEBUG("(%d:%d:%s) calling multiply on index %d (%d:%d:%d)\n",
              tier, index, bitIndex.to_string().c_str(), (index << 3) | productIndex, i, j, k);
          child[childIndex]->matmul(*AInfo->child[childIndexA],
              *BInfo->child[childIndexB],
              (index << 3) | productIndex,
              CkCallback(CkIndex_Node::matmulDone(indexMsg), thisProxy));
        }
      }
    }
  }

  LOG_DEBUG("(%d:%d:%s) done\n", tier, index, bitIndex.to_string().c_str());
}

void Node::matmulDone (IntMsg *index)
{
  std::bitset<20> bitIndex = std::bitset<20>(matmulIndex);
  int kIndex = (index->i >> 3) & 1;
  int CIndex = (index->i & 6) >> 1;
  int tierIndex = index->i & 7;
  LOG_DEBUG("(%d:%d:%s) matmulDone called with index = %d, CIndex = %d, "
      "kIndex = %d, tierIndex = %d\n",
      tier, matmulIndex, bitIndex.to_string().c_str(), index->i, CIndex,
      kIndex, tierIndex);
  for(int i = 0; i < numberChildProducts; i++)
  {
    if(childProducts[i] == index->i)
    {
      for(int j = i+1; j < numberChildProducts; j++)
      {
        childProducts[j-1] = childProducts[j];
      }
      numberChildProducts--;
      break;
    }
  }

  if(numberChildProducts == 0)
  {
    LOG_DEBUG("(%d:%d:%s) matmulDone sending to parent\n", tier, matmulIndex,
        bitIndex.to_string().c_str());
    parentDone.send(new IntMsg(matmulIndex));
  }
}

#include "Node.def.h"
