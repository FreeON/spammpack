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
  this->tier = tier;
}

/** Get information on node.
 *
 * @return A NodeMsg object.
 */
NodeMsg * Node::info ()
{
  return new NodeMsg (iLower, iUpper, jLower, jUpper, blocksize, child, data);
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
      *child[childIndex] = CProxy_Node::ckNew(tier+1, blocksize,
          newILower, newJLower,
          newILower+width, newJLower+width);
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
#ifdef CALLBACK
void Node::matmul (CProxy_Node A, CProxy_Node B, int productIndex, CkCallback &done)
#elif defined(FUTURES)
void Node::matmul (CProxy_Node A, CProxy_Node B, int productIndex, CkFuture f)
#else
#error "FIXME"
#endif
{
#ifdef CALLBACK
  /* Store callback. */
  parentDone = done;

  /* Construct an integer message. */
  IntMsg *indexMsg = new IntMsg(productIndex);
#endif
#ifdef FUTURES
  IntMsg *result = new IntMsg();
#endif

  /* The width of the C matrix block. */
  int width = iUpper-iLower;

  int counter = 0;

  NodeMsg *AInfo = A.info();
  NodeMsg *BInfo = B.info();

  std::string tagstr = getTagString(tier, productIndex);

  if(width == blocksize)
  {
    if(AInfo->data != NULL && BInfo->data != NULL)
    {
      if(data == NULL)
      {
        data = new double[blocksize*blocksize];
        memset(data, 0, blocksize*blocksize*sizeof(double));
      }

      /* Get matrix blocks from A and B. */
      DataMsg *AData = A.getData();
      DataMsg *BData = B.getData();

      if(iLower < 0 || jLower < 0 || iUpper < 0 || jUpper < 0)
      {
        LOG_INFO("negative\n");
        CkExit();
      }

//#define BLOCK_MULTIPLY

#ifdef BLOCK_MULTIPLY
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
#else

#define BLOCK_FAKE_SKIP
//#define BLOCK_FAKE_AB_READ
//#define BLOCK_FAKE_AB_READ_MULTIPLY
//#define BLOCK_FAKE_C_READ
//#define BLOCK_FAKE_C_READ_WRITE

#ifndef BLOCK_FAKE_SKIP
      for(int i = iLower; i < iUpper; i++) {
        for(int j = jLower; j < jUpper; j++) {
          for(int k = AInfo->jLower; k < AInfo->jUpper; k++)
          {
#ifdef BLOCK_FAKE_AB_READ
            /* Read matrix elements of A and B. */
            if(AData->data[blockIndex(i, k, AInfo->iLower, AInfo->jLower, AInfo->blocksize)] < 0) { printf("A\n"); }
            if(BData->data[blockIndex(k, j, BInfo->iLower, BInfo->jLower, BInfo->blocksize)] < 0) { printf("B\n"); }
#elif defined(BLOCK_FAKE_AB_READ_MULTIPLY)
            /* Read and multiply matrix elements of A and B. */
            if(AData->data[blockIndex(i, k, AInfo->iLower, AInfo->jLower, AInfo->blocksize)]
                *BData->data[blockIndex(k, j, BInfo->iLower, BInfo->jLower, BInfo->blocksize)] < 0)
            {
              printf("B\n");
            }
#endif
#ifdef BLOCK_FAKE_C_READ
            /* Read matrix elements of C. */
            if(data[blockIndex(i, j, iLower, jLower, blocksize)] < 0) { printf("C\n"); }
#elif defined(BLOCK_FAKE_C_READ_WRITE)
            /* Read and write matrix elements of C. */
            data[blockIndex(i, j, iLower, jLower, blocksize)]++;
#endif
          }
        }
      }
#endif

#endif

      delete AData;
      delete BData;

      result->i = 1;
    }

    else
    {
      LOG_INFO("%s [FIXME] skipping C(%d:%d,%d:%d) += A(%d:%d,%d:%d)*B(%d:%d,%d:%d), delete C\n",
          iLower, iUpper, jLower, jUpper,
          AInfo->iLower, AInfo->iUpper, AInfo->jLower, AInfo->jUpper,
          BInfo->iLower, BInfo->iUpper, BInfo->jLower, BInfo->jUpper,
          tagstr.c_str());
    }
#ifdef CALLBACK
    /* Signal that we are done. */
    done.send(indexMsg);
#endif
  }

  else
  {

#if defined(FUTURES)
    CkFuture product_f[8];
    int numberProducts = 0;
#endif
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

#ifdef CALLBACK
          int kIndex = productIndex & 1;

          queuedProducts[kIndex][numberQueued[kIndex]++] = childProductIndex;

          LOG_DEBUG("%s added %d to queue[%d], has %d elements now\n",
              tagstr.c_str(), childProductIndex, kIndex, numberQueued[kIndex]);
#endif
          if(AInfo->child[childIndexA] == NULL || BInfo->child[childIndexB] == NULL)
          {
            LOG_INFO("%s [FIXME] delete C (product %d is NULL).\n",
                tagstr.c_str(), childProductIndex);
#ifdef CALLBACK
            thisProxy.matmulDone(new IntMsg(childProductIndex));
#endif
            continue;
          }

          if(child[childIndex] == NULL)
          {
            child[childIndex] = new CProxy_Node;
            *child[childIndex] = CProxy_Node::ckNew(tier+1, blocksize,
                iLower+width/2*i, jLower+width/2*j,
                iLower+width/2*(i+1), jLower+width/2*(j+1));
          }

          Counter::increment();
#if defined(CALLBACK)
          child[childIndex]->matmul(*AInfo->child[childIndexA],
              *BInfo->child[childIndexB],
              childProductIndex,
              CkCallback(CkIndex_Node::matmulDone(indexMsg), thisProxy));
#elif defined(FUTURES)
          product_f[numberProducts] = CkCreateFuture();
          child[childIndex]->matmul(*AInfo->child[childIndexA],
              *BInfo->child[childIndexB],
              childProductIndex,
              product_f[numberProducts]);
          numberProducts++;
          if(numberProducts > 8)
          {
            LOG_ERROR("numberProducts = %d\n", numberProducts);
            CkExit();
          }
#else
#error "FIXME"
#endif
        }
      }
    }
#ifdef FUTURES
    counter++;
    for(int i = 0; i < numberProducts; i++)
    {
      IntMsg *m = (IntMsg*) CkWaitFuture(product_f[i]);
      result->i += m->i;
      delete m;
      CkReleaseFuture(product_f[i]);
    }
#endif
  }

  delete AInfo;
  delete BInfo;

#if defined(FUTURES)
  CkSendToFuture(f, result);
#endif
}

#ifdef CALLBACK
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
#endif

#include "Node.def.h"
