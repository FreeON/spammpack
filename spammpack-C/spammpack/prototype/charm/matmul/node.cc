/** @file
 *
 * The implementation of the Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "node.h"
#include "messages.h"
#include "logger.h"
#include "index.h"

#include <bitset>

/** The constructor.
 *
 * @param N The matrix size.
 * @param depth The matrix depth.
 * @param blocksize The blocksize.
 * @param tier The tier this node is on.
 * @param seed The random number seed.
 */
Node::Node (int N, int depth, int blocksize, int tier, unsigned int seed)
{
  this->N = N;
  this->depth = depth;
  this->blocksize = blocksize;
  this->tier = tier;
  this->seed = seed;

  this->iLower = thisIndex.x*blocksize;
  this->iUpper = (thisIndex.x+1)*blocksize;
  this->jLower = thisIndex.y*blocksize;
  this->jUpper = (thisIndex.y+1)*blocksize;

  this->block = new double[blocksize*blocksize];
  memset(this->block, 0, sizeof(double)*blocksize*blocksize);

  for(int i = iLower; i < iUpper && i < N; i++) {
    for(int j = jLower; j < jUpper && j < N; j++)
    {
      this->block[BLOCK_INDEX(i, j, iLower, jLower, blocksize)] =
        rand_r(&seed)/(double) RAND_MAX;
    }
  }

#ifdef DEBUG_OUTPUT
  /* For debugging. */
  printDense(blocksize, this->block, "tier %d, Node(%d,%d), initial matrix",
      tier, thisIndex.x, thisIndex.y);
#endif

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

  DEBUG("tier %d, Node(%d,%d), index %s, constructing\n", tier, thisIndex.x,
      thisIndex.y, toBinary(index).c_str());
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
Node::Node (CkMigrateMessage *msg)
{
  DEBUG("tier %d, Node(%d,%d) migration constructor\n", -1, thisIndex.x, thisIndex.y);
}

/** The destructor.
 */
Node::~Node (void)
{
  DEBUG("tier %d, Node(%d,%d) destructor\n", tier, thisIndex.x, thisIndex.y);
  delete[] block;
}

/** The PUP method.
 */
void Node::pup (PUP::er &p)
{
  CBase_Node::pup(p);

  p|seed;
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

  int numberElements = (block == NULL ? 0 : blocksize*blocksize);
  p|numberElements;

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

  DEBUG("tier %d, Node(%d,%d) pup()\n", tier, thisIndex.x, thisIndex.y);
}

/** Get information on a Node.
 *
 * @return A message.
 */
NodeInfoMsg * Node::info (void)
{
  DEBUG("Node(%d,%d) getting node info on index %s\n",
      thisIndex.x, thisIndex.y, toBinary(index).c_str());

  return new NodeInfoMsg(index, norm, norm_2);
}

/** Get the dense submatrix block.
 *
 * @return The submatrix block.
 */
DenseMatrixMsg * Node::getBlock (void)
{
  DEBUG("Node(%d,%d) getting block\n", thisIndex.x, thisIndex.y);
  DenseMatrixMsg *m = new (blocksize*blocksize) DenseMatrixMsg();
  if(block != NULL)
  {
    memcpy(m->A, block, sizeof(double)*blocksize*blocksize);
  }
  else
  {
    memset(m->A, 0, sizeof(double)*blocksize*blocksize);
  }
  return m;
}

/** Add a submatrix block to this Node.
 *
 * @param blocksize The blocksize.
 * @param A The dense matrix.
 */
void Node::add (int blocksize, double *A)
{
  DEBUG("tier %d, Node(%d,%d) Adding back to C\n", tier, thisIndex.x, thisIndex.y);
  for(int i = 0; i < blocksize*blocksize; i++)
  {
    block[i] += A[i];
  }

#ifdef DEBUG_OUTPUT
  /* For debugging. */
  printDense(blocksize, block, "tier %d, Node(%d,%d) Adding back to C", tier,
      thisIndex.x, thisIndex.y);
#endif
}

#include "node.def.h"
