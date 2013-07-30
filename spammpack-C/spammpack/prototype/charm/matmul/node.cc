/** @file
 *
 * The implementation of the Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "node.h"
#include "logger.h"
#include "messages.h"
#include "index.h"
#include <bitset>
#include <string>

/** The constructor.
 */
Node::Node ()
{
  block = NULL;
}

/** The constructor.
 */
Node::Node (int N, int depth, int blocksize, int tier)
{
  this->N = N;
  this->depth = depth;
  this->blocksize = blocksize;
  this->tier = tier;

  this->iLower = thisIndex.x*blocksize;
  this->iUpper = (thisIndex.x+1)*blocksize;
  this->jLower = thisIndex.y*blocksize;
  this->jUpper = (thisIndex.y+1)*blocksize;

  this->tierNodeSet = false;

  this->block = NULL;

  /* Calculate the linear index. */
  std::bitset<8*sizeof(int)> iIndex(thisIndex.x);
  std::bitset<8*sizeof(int)> jIndex(thisIndex.y);
  std::bitset<8*sizeof(int)> tempIndex(1);
  for(int i = 0; i < tier; i++)
  {
    tempIndex <<= 2;
    if(iIndex[i]) { tempIndex |= 2; }
    if(jIndex[i]) { tempIndex |= 1; }
  }
  index = tempIndex.to_ulong();

  DEBUG("tier %d, index %d, constructing\n", tier, index);
  DEBUG("Node(%d,%d) constructor\n", thisIndex.x, thisIndex.y);
}

/** The migration constructor. */
Node::Node (CkMigrateMessage *msg)
{
  DEBUG("Node(%d,%d) migration constructor\n", thisIndex.x, thisIndex.y);
}

/** The destructor.
 */
Node::~Node ()
{
  DEBUG("Node(%d,%d) destructor\n", thisIndex.x, thisIndex.y);
  delete[] block;
}

/** The PUP method.
 */
void Node::pup (PUP::er &p)
{
  CBase_Node::pup(p);
  p|N;
  p|blocksize;
  p|depth;
  p|tier;
  p|iLower;
  p|iUpper;
  p|jLower;
  p|jUpper;
  p|index;
  p|norm;
  p|norm_2;
  p|tierNodeSet;

  if(tierNodeSet)
  {
    p|tierNode;
  }

  int numberElements = (block == NULL ? 0 : blocksize*blocksize);
  p|numberElements;

  if(p.isUnpacking())
  {
    DEBUG("pup: Node(%d,%d) unpacking %d elements\n", thisIndex.x, thisIndex.y, numberElements);
  }
  else
  {
    if(p.isSizing())
    {
      DEBUG("pup: Node(%d,%d) sizing %d elements\n", thisIndex.x, thisIndex.y, numberElements);
    }
    else
    {
      DEBUG("pup: Node(%d,%d) packing %d elements\n", thisIndex.x, thisIndex.y, numberElements);
    }
  }

  if(numberElements > 0)
  {
    if(p.isUnpacking())
    {
      block = new double[numberElements];
    }
    PUParray(p, block, numberElements);
  }
  else
  {
    if(p.isUnpacking()) { block = NULL; }
  }
}

/** Get the dense submatrix block.
 *
 * @return The submatrix block.
 */
NodeBlockMsg * Node::getBlock ()
{
  DEBUG("Node(%d,%d) getting block\n", thisIndex.x, thisIndex.y);
  NodeBlockMsg *m = new (blocksize*blocksize) NodeBlockMsg();
  memcpy(m->block, block, sizeof(double)*blocksize*blocksize);
  return m;
}

/** Get a matrix element.
 *
 * @param i The row index.
 * @param j The column index.
 *
 * @return The matrix element.
 */
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
      if(i < iLower || j < jLower || i >= iUpper || j >= jUpper)
      {
        ABORT("out of bounds\n");
      }
      DEBUG("found block[%d] = %f (block norm = %f)\n",
          BLOCK_INDEX(i, j, iLower, jLower, blocksize),
          block[BLOCK_INDEX(i, j, iLower, jLower, blocksize)],
          norm);
      return new DoubleMsg(block[BLOCK_INDEX(i, j, iLower, jLower, blocksize)]);
    }
  }

  else
  {
    ABORT("not at lowest tier\n");
  }
}

/** Get information on a Node.
 *
 * @return A message.
 */
NodeInfoMsg * Node::info ()
{
  DEBUG("getting node info on index %d\n", index);
  NodeInfoMsg *msg = new NodeInfoMsg(index, norm, norm_2);
  return msg;
}

/** Print the PE this Node is on. */
void Node::printPE (CkCallback &cb)
{
  CkPrintf("node is on PE %d\n", CkMyPe());
  contribute(cb);
}

/** Initialize a Matrix with random numbers. */
void Node::random (CkCallback &cb)
{
  DEBUG("generating random matrix\n");
  thisProxy.initialize(initRandom, cb);
}

/** Initialize a Matrix with zeros. */
void Node::zero (CkCallback &cb)
{
  DEBUG("setting matrix to zero\n");
  thisProxy.initialize(initZero, cb);
}

/** Set a matrix block directly.
 *
 * @param numberElements The number of elements contained in A.
 * @param A The dense submatrix.
 * @param cb The callback.
 */
void Node::set (int numberElements, double *A, CkCallback &cb)
{
  if(numberElements != blocksize*blocksize)
  {
    ABORT("blocksize mismatch\n");
  }

  if(block == NULL)
  {
    block = new double[blocksize*blocksize];
  }

  memcpy(block, A, numberElements*sizeof(double));

#ifdef DEBUG_OUTPUT
  DEBUG("Node(%d,%d) setting block\n", thisIndex.x, thisIndex.y);
  for(int i = 0; i < blocksize; i++) {
    for(int j = 0; j < blocksize; j++)
    {
      CkPrintf(" %e", block[BLOCK_INDEX(i, j, 0, 0, blocksize)]);
    }
    CkPrintf("\n");
  }
#endif

  norm_2 = 0;
  for(int i = 0; i < numberElements; i++)
  {
    norm_2 += block[i]*block[i];
  }
  norm = sqrt(norm_2);

  cb.send();
}

/** Initialize a Node.
 *
 * @param initType How to initialize the Matrix.
 */
void Node::initialize (int initType, CkCallback &cb)
{
  if(tier < depth)
  {
    ABORT("not at depth\n");
  }

  DEBUG("Node(%d,%d) index %d, initializing\n", thisIndex.x, thisIndex.y, index);

  if(block == NULL)
  {
    block = new double[blocksize*blocksize];
    memset(block, 0, sizeof(double)*blocksize*blocksize);
  }

  switch(initType)
  {
    case initRandom:
      norm_2 = 0;
      for(int i = 0; i < blocksize && i+iLower < N; i++) {
        for(int j = 0; j < blocksize && j+jLower < N; j++)
        {
          double aij = rand()/(double) RAND_MAX;
          block[BLOCK_INDEX(i, j, 0, 0, blocksize)] = aij;
          norm_2 += aij*aij;
        }
      }
      norm = sqrt(norm_2);
      break;

    case initZero:
      norm = 0;
      norm_2 = 0;
      break;

    default:
      ABORT("unknown initType\n");
      break;
  }

  cb.send();
}

/** Print the PEs the leafs sit on. */
void Node::printLeafPes (CkCallback &cb)
{
  std::string bitString = std::bitset<8*sizeof(unsigned int)>(index).to_string();
  while(bitString[0] == '0')
  {
    bitString.erase(0, 1);
  }
  CkPrintf("leaf %s on PE %d\n", bitString.c_str(), CkMyPe());
  contribute(cb);
}

/** Add a submatrix block to this Node.
 *
 * @param A The dense matrix.
 */
void Node::add (int blocksize, double *A)
{
  DEBUG("Adding block to A(%d,%d)\n", thisIndex.x, thisIndex.y);

  for(int i = 0; i < blocksize*blocksize; i++)
  {
    block[i] += A[i];
  }
}

/** Update norms of a tier. */
void Node::updateNorms (CkCallback &cb)
{
  /* Get norms from nodes underneath this one. */
  norm_2 = 0;
  for(int i = 0; i < 2; i++) {
    for(int j = 0; j < 2; j++)
    {
      NodeInfoMsg *childInfo = tierNode((thisIndex.x << 1)+i, (thisIndex.y << 1)+j).info();
      norm_2 += childInfo->norm_2;
      delete childInfo;
    }
  }
  norm = sqrt(norm_2);
  DEBUG("tier %d, Node(%d,%d) setting norm = %f\n", tier, thisIndex.x, thisIndex.y, norm);
  contribute(cb);
}

/** Set the tierNode
 *
 * @param tierNode The node list for the next tier.
 */
void Node::setTierNode (CProxy_Node tierNode, CkCallback &cb)
{
  this->tierNodeSet = true;
  this->tierNode = tierNode;
  contribute(cb);
}

#include "node.def.h"
