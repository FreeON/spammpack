#include "Node.h"
#include "Utilities.h"

#include <bitset>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

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
  this->numberQueued[0] = 0;
  this->numberQueued[1] = 0;
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

std::string getTagString (int tier, int productIndex)
{
  std::bitset<20> bitIndex = std::bitset<20>(productIndex);
  std::ostringstream tag;
  tag << "(" << tier << ":" << productIndex << ":" << bitIndex.to_string() << ")";
  return tag.str();
}

/** Multiply two matrices.
 *
 * @param A Node A.
 * @param B Node B.
 * @param productIndex The 3D linear index of the product contribution.
 * @param done The callback of the parent product contribution.
 */
void Node::matmul (CProxy_Node A, CProxy_Node B, int productIndex, CkCallback &done)
{
  /* Store callback. */
  parentDone = done;

  /* Construct an integer message. */
  IntMsg *indexMsg = new IntMsg(productIndex);

  /* The width of the C matrix block. */
  int width = iUpper-iLower;

  NodeMsg *AInfo = A.info();
  NodeMsg *BInfo = B.info();

  std::string tagstr = getTagString(tier, productIndex);

  LOG_DEBUG("%s starting multiply, C(%d:%d,%d:%d) += A(%d:%d,%d:%d)*B(%d:%d,%d:%d)\n",
      tagstr.c_str(),
      iLower+1, iUpper, jLower+1, jUpper,
      AInfo->iLower+1, AInfo->iUpper, AInfo->jLower+1, AInfo->jUpper,
      BInfo->iLower+1, BInfo->iUpper, BInfo->jLower+1, BInfo->jUpper);

  if(width == blocksize)
  {
    DataMsg *AData = A.getData();
    DataMsg *BData = B.getData();

    if(AData->data != NULL && BData->data != NULL)
    {
      if(data == NULL)
      {
        data = new double[blocksize*blocksize];
        memset(data, 0, blocksize*blocksize*sizeof(double));
      }

      LOG_DEBUG("%s block multiply\n", tagstr.c_str());
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
    }

    else
    {
      LOG_ERROR("%s [FIXME] delete C\n", tagstr.c_str());
    }

    /* Signal that we are done. */
    LOG_DEBUG("%s sending %d to done\n", tagstr.c_str(), productIndex);
    done.send(indexMsg);
  }

  else
  {
    LOG_DEBUG("%s descending...\n", tagstr.c_str());

    for(int i = 0; i < 2; i++) {
      for(int j = 0; j < 2; j++)
      {
        int childIndex = (i << 1) | j;
        for(int k = 0; k < 2; k++)
        {
          int childIndexA = (i << 1) | k;
          int childIndexB = (k << 1) | j;

          int convolutionIndex = (i << 2) | (j << 1) | k;
          int childProductIndex = (productIndex << 3) | convolutionIndex;


          int kIndex = productIndex & 1;
          queuedProducts[kIndex][numberQueued[kIndex]++] = childProductIndex;

          LOG_DEBUG("%s added %d to queue[%d], has %d elements now\n",
              tagstr.c_str(), childProductIndex, kIndex, numberQueued[kIndex]);

          if(AInfo->child[childIndexA] == NULL || BInfo->child[childIndexB] == NULL)
          {
            LOG_ERROR("%s [FIXME] delete C (product %d is NULL).\n",
                tagstr.c_str(), childProductIndex);
            thisProxy.matmulDone(new IntMsg(childProductIndex));
            continue;
          }
          if(child[childIndex] == NULL)
          {
            child[childIndex] = new CProxy_Node;
            *child[childIndex] = CProxy_Node::ckNew(tier+1, blocksize,
                iLower+width/2*i, jLower+width/2*j, iLower+width/2*(i+1),
                jLower+width/2*(j+1));
          }
          LOG_DEBUG("%s calling multiply on index %d (%d:%d:%d)\n",
              tagstr.c_str(), childProductIndex, i, j, k);
          child[childIndex]->matmul(*AInfo->child[childIndexA],
              *BInfo->child[childIndexB],
              childProductIndex,
              CkCallback(CkIndex_Node::matmulDone(indexMsg), thisProxy));
        }
      }
    }
  }

  LOG_DEBUG("%s done\n", tagstr.c_str());
}

void Node::matmulDone (IntMsg *productIndex)
{
  std::string tagstr = getTagString(tier, productIndex->i >> 3);

  int kIndex = (productIndex->i >> 3) & 1;
  int CIndex = (productIndex->i & 6) >> 1;
  int tierIndex = productIndex->i & 7;

  LOG_DEBUG("%s matmulDone called with index = %d, CIndex = %d, "
      "kIndex = %d, tierIndex = %d, numberQueued = %i\n",
      tagstr.c_str(), productIndex->i, CIndex, kIndex, tierIndex,
      numberQueued[(productIndex->i >> 3) & 1]);

  for(int i = 0; i < numberQueued[kIndex]; i++)
  {
    if(queuedProducts[kIndex][i] == productIndex->i)
    {
      for(int j = i+1; j < numberQueued[kIndex]; j++)
      {
        queuedProducts[kIndex][j-1] = queuedProducts[kIndex][j];
      }
      numberQueued[kIndex]--;
      LOG_DEBUG("%s removed %d from queue\n", tagstr.c_str(), productIndex->i);
      break;
    }
  }

  if(numberQueued[kIndex] == 0)
  {
    LOG_DEBUG("%s matmulDone sending %d to parent\n", tagstr.c_str(), productIndex->i >> 3);
    parentDone.send(new IntMsg(productIndex->i >> 3));
  }
}

#include "Node.def.h"
