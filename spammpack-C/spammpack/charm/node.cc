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
#include "types.h"
#include "utilities.h"

#include <assert.h>
#include <bitset>

/** Some convenience macros for logging. Wrap the logging format string with
 * LB and LE. */
#define LB "tier %d Node(%d,%d) "

/** Some convenience macros for logging. Wrap the logging format string with
 * LB and LE. */
#define LE , tier, thisIndex.x, thisIndex.y

/** The constructor.
 *
 * @param N The matrix size.
 * @param depth The matrix depth.
 * @param blocksize The blocksize.
 * @param tier The tier this node is on.
 */
Node::Node (int N, int depth, int blocksize, int tier)
{
  this->N = N;
  this->blocksize = blocksize;
  this->depth = depth;
  this->tier = tier;

  this->iLower = thisIndex.x*blocksize;
  this->iUpper = (thisIndex.x+1)*blocksize;
  this->jLower = thisIndex.y*blocksize;
  this->jUpper = (thisIndex.y+1)*blocksize;

  this->norm = 0;
  this->norm_2 = 0;

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

  DEBUG(LB"index %s, constructing\n"LE, toBinary(index).c_str());
}

/** The migration constructor.
 *
 * @param msg The migration message.
 */
Node::Node (CkMigrateMessage *msg)
{
  DEBUG("Node(%d,%d) migration constructor\n", thisIndex.x, thisIndex.y);
}

/** The destructor.
 */
Node::~Node (void)
{
  DEBUG(LB"destructor\n"LE);
  delete[] block;
}

/** The PUP method.
 *
 * @param p The Pup::er object.
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

  DEBUG(LB"pup()\n"LE);
}

/** Get information on a Node.
 *
 * @return A message.
 */
NodeInfoMsg * Node::info (void)
{
  DEBUG(LB"getting node info on index %s, "
      "norm = %e, norm_2 = %e\n"LE,
      toBinary(index).c_str(), norm, norm_2);

  return new NodeInfoMsg(index, norm, norm_2);
}

/** Get the dense submatrix block.
 *
 * @return The submatrix block.
 */
DenseMatrixMsg * Node::getBlock (void)
{
  DEBUG(LB"getting block\n"LE);

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

/** Set a matrix block in this Node.
 *
 * @param blocksize The blocksize.
 * @param A The matrix.
 * @param cb The callback to send to.
 */
void Node::set (int blocksize, double *A, CkCallback &cb)
{
  assert(blocksize == this->blocksize);

  if(block == NULL)
  {
    block = new double[blocksize*blocksize];
  }
  memcpy(block, A, sizeof(double)*blocksize*blocksize);

  norm_2 = 0;
  for(int i = 0; i < blocksize*blocksize; i++)
  {
    norm_2 += block[i]*block[i];
  }
  norm = sqrt(norm_2);

  DEBUG(LB"norm = %e\n"LE, norm);

#ifdef DEBUG_OUTPUT
  printDense(blocksize, block, LB"setting block:"LE);
#endif

  cb.send();
}

/** Set the norm of this Node based on the norms of the Node objects
 * underneath it.
 *
 * @param nodes The Node chare array underneath this tier.
 * @param cb The callback to reduce to.
 */
void Node::setNorm (CProxy_Node nodes, CkCallback &cb)
{
  DEBUG(LB"updating norms\n"LE);

  norm_2 = 0;
  for(int i = 0; i < 2; i++)
  {
    int nextX = (thisIndex.x << 1) | i;
    for(int j = 0; j < 2; j++)
    {
      int nextY = (thisIndex.y << 1) | j;
      NodeInfoMsg *msg = nodes(nextX, nextY).info();
      DEBUG(LB"got tier %d, Node(%d,%d) norm = %e\n"LE, tier+1, nextX, nextY,
          norm);
      norm_2 += msg->norm_2;
      delete msg;
    }
  }
  norm = sqrt(norm_2);
  DEBUG(LB"norm = %e\n"LE, norm);
  contribute(cb);
}

/** Add a submatrix block to this Node.
 *
 * @param blocksize The blocksize.
 * @param A The dense matrix.
 */
void Node::blockAdd (int blocksize, double *A)
{
  assert(blocksize == this->blocksize);

  DEBUG(LB"Adding back to C\n"LE);
  if(block == NULL)
  {
    block = new double[blocksize*blocksize];
    memset(block, 0, sizeof(double)*blocksize*blocksize);
  }

  norm_2 = 0;
  for(int i = 0; i < blocksize*blocksize; i++)
  {
    block[i] += A[i];
    norm_2 += block[i]*block[i];
  }
  norm = sqrt(norm_2);

#ifdef DEBUG_OUTPUT
  /* For debugging. */
  printDense(blocksize, block, LB"Adding back to C"LE);
#endif
}

/** Add to this matrix a second matrix scaled by a factor, i.e.
 *
 * @f[ A \leftarrow \alpha A + \beta B @f]
 *
 * @param alpha Factor @f$ \alpha @f$.
 * @param beta Factor @f$ \beta @f$.
 * @param B Matrix B.
 * @param cb The callback for the reduction.
 */
void Node::add (double alpha, double beta, CProxy_Node B, CkCallback &cb)
{
  DEBUG(LB"adding, alpha = %e\n"LE, alpha);
  DenseMatrixMsg *BData = B(thisIndex.x, thisIndex.y).getBlock();
  if(block == NULL)
  {
    block = new double[blocksize*blocksize];
    memset(block, 0, sizeof(double)*blocksize*blocksize);
  }
  for(int i = 0; i < blocksize*blocksize; i++)
  {
    block[i] = alpha*block[i]+beta*BData->A[i];
  }
  delete BData;
  contribute(cb);
}

/** Reduce the trace of the matrix.
 *
 * @param cb The reduction target.
 */
void Node::trace (CkCallback &cb)
{
  double trace = 0;

  if(thisIndex.x == thisIndex.y)
  {
    for(int i = 0; i < blocksize; i++)
    {
      trace += block[BLOCK_INDEX(i, i, 0, 0, blocksize)];
    }
  }
  contribute(sizeof(double), &trace, CkReduction::sum_double, cb);
}

/** Create a PE map of the Matrix @link Node nodes @endlink.
 *
 * @param cb The callback for the reduction.
 */
void Node::PEMap (CkCallback &cb)
{
  DEBUG(LB"PE %d\n"LE, CkMyPe());

  struct PEMap_Node_t *result = (struct PEMap_Node_t*) malloc(sizeof(struct PEMap_Node_t));

  result->index[0] = thisIndex.x;
  result->index[1] = thisIndex.y;
  result->PE = CkMyPe();
  result->norm = norm;

  contribute(sizeof(struct PEMap_Node_t), result, CkReduction::set, cb);
}

#include "node.def.h"