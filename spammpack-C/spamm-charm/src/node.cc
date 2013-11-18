/** @file
 *
 * The implementation of the Node class.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"

#include "node.h"

#include "chunk.h"
#include "index.h"
#include "logger.h"
#include "messages.h"
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
Node::Node (int N, int depth, int blocksize, int N_basic, int tier)
{
  this->N = N;
  this->blocksize = blocksize;
  this->N_basic = N_basic;
  this->depth = depth;
  this->tier = tier;

  this->iLower = thisIndex.x*blocksize;
  this->iUpper = (thisIndex.x+1)*blocksize;
  this->jLower = thisIndex.y*blocksize;
  this->jUpper = (thisIndex.y+1)*blocksize;

  this->norm = 0;
  this->norm_2 = 0;

  this->chunksize = 0;
  this->chunk = NULL;

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
  chunk = NULL;
}

/** The destructor.
 */
Node::~Node (void)
{
  DEBUG(LB"destructor\n"LE);
  if(chunk != NULL)
  {
    DEBUG("free'ing chunk at %p\n", chunk);
    free(chunk);
    chunk = NULL;
  }
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
  p|N_basic;
  p|depth;
  p|tier;
  p|iLower;
  p|iUpper;
  p|jLower;
  p|jUpper;
  p|index;
  p|norm;
  p|norm_2;
  p|chunksize;
  bool chunkSet = (chunk != NULL);
  p|chunkSet;
  if(chunkSet)
  {
    if(p.isUnpacking())
    {
      chunk = malloc(chunk_sizeof(blocksize, N_basic));
      assert(chunksize == chunk_sizeof(blocksize, N_basic));
    }
    PUParray(p, (char*) chunk, chunksize);
  }
  DEBUG(LB"pup()\n"LE);
}

/** Intialize the Node.
 *
 * This is merely a convenience method to force the Charm++ runtime to
 * instantiate the Node. The RTS is sometimes a little too lazy...
 *
 * @param cb The callback to send to when done.
 */
void Node::init (CkCallback &cb)
{
  INFO(LB"initializing\n"LE);
  contribute(cb);
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

  return new NodeInfoMsg(index, iLower, iUpper, jLower, jUpper, norm, norm_2);
}

/** Get the dense submatrix block.
 *
 * @return The submatrix block.
 */
DenseMatrixMsg * Node::toDense (void)
{
  DEBUG(LB"getting block\n"LE);

  DenseMatrixMsg *msg = new (blocksize*blocksize) DenseMatrixMsg(blocksize, blocksize);

  if(chunk != NULL)
  {
    double *A = chunk_to_dense(chunk);
    memcpy(msg->A, A, sizeof(double)*blocksize*blocksize);
    free(A);
  }

  return msg;
}

/** Get the matrix Chunk.
 *
 * @return The matrix Chunk.
 */
ChunkMsg * Node::getChunk (void)
{
  DEBUG(LB"chunk at %p, chunksize = %d\n"LE, chunk, chunksize);
  ChunkMsg *msg;
  if(chunk == NULL)
  {
    msg = new (0) ChunkMsg(chunksize, (char*) chunk);
  }

  else
  {
    msg = new (chunksize) ChunkMsg(chunksize, (char*) chunk);
  }

  return msg;
}

/** Set a matrix block in this Node.
 *
 * @param blocksize The blocksize.
 * @param A The matrix.
 * @param cb The callback to send to.
 */
void Node::set (int blocksize, double *A, CkCallback &cb)
{
  assert(tier == depth);
  assert(blocksize == this->blocksize);

  if(chunk == NULL)
  {
    chunk = chunk_alloc(blocksize, N_basic, N, iLower, jLower);
    chunksize = chunk_sizeof(blocksize, N_basic);
    DEBUG(LB"creating new Chunk at %p, chunksize = %llu\n"LE, chunk, chunksize);
  }
  chunk_set(chunk, A);
  norm_2 = chunk_get_norm(chunk);
  norm = sqrt(norm_2);

#ifdef DEBUG_OUTPUT
  chunk_print(chunk, LB"setting block:"LE);
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
  assert(tier < depth);

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

/** Add a scaled submatrix block to this Node.
 *
 * @f[ A \leftarrow \alpha A @f]
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param blocksize The blocksize.
 * @param A The dense matrix.
 */
void Node::chunkAdd (double alpha, size_t chunksize, char *chunk)
{
  assert(tier == depth);
  assert(blocksize == this->blocksize);
  assert(chunksize == chunk_sizeof(blocksize, N_basic));

  if(this->chunk == NULL)
  {
    this->chunk = chunk_alloc(blocksize, N_basic, N, iLower, jLower);
    this->chunksize = chunksize;
    DEBUG(LB"creating new Chunk at %p\n"LE, this->chunk);
  }

  DEBUG(LB"Adding back to C with Chunk at %p and alpha = %e\n"LE, this->chunk, alpha);
  chunk_add(1, this->chunk, alpha, chunk);

#ifdef DEBUG_OUTPUT
  /* For debugging. */
  chunk_print(chunk, LB"Adding back to C"LE);
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
  assert(tier == depth);

  DEBUG(LB"alpha = %e, beta = %e, chunk at %p\n"LE, alpha, beta, chunk);
  ChunkMsg *BChunk = B(thisIndex.x, thisIndex.y).getChunk();

  if(chunk == NULL)
  {
    DEBUG(LB"allocating new chunk\n"LE);
    chunk = chunk_alloc(blocksize, N_basic, N, iLower, jLower);
  }

  if(BChunk->chunksize != 0)
  {
    chunk_add(alpha, chunk, beta, BChunk->chunk);
  }

  else
  {
    chunk_scale(alpha, chunk);
  }

  delete BChunk;

  contribute(cb);
}

/** Reduce the trace of the matrix.
 *
 * @param cb The reduction target.
 */
void Node::trace (CkCallback &cb)
{
  double trace = 0;

  assert(tier == depth);

  if(iLower == jLower && iUpper == jUpper)
  {
    if(chunk != NULL)
    {
      trace = chunk_trace(chunk);
    }
  }
  DEBUG(LB"trace(%d,%d) = %e\n"LE, thisIndex.x, thisIndex.y, trace);

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

/** Scale a Node by a scalar factor.
 *
 * @f[ A \leftarrow \alpha A @f].
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param cb The reduction target callback.
 */
void Node::scale (double alpha, CkCallback &cb)
{
  assert(tier == depth);

  if(chunk != NULL)
  {
    DEBUG(LB"alpha = %e of Chunk at %p\n"LE, alpha, chunk);
    chunk_scale(alpha, chunk);
    norm_2 *= alpha*alpha;
    norm = sqrt(norm_2);
  }

  else
  {
    DEBUG(LB"block == NULL\n"LE);
  }

  contribute(cb);
}
/** Add a the scaled identity matrix. This operation affects only those Nodes
 * which are on the diagonal.
 *
 * @f[ A \leftarrow A + \alpha I @f]
 *
 * @param alpha The scalar alpha.
 * @param cb The reduction target callback.
 */
void Node::addIdentity (double alpha, CkCallback &cb)
{
  assert(tier == depth);

  DEBUG(LB"alpha = %e\n"LE, alpha);

  if(chunk != NULL)
  {
    if(iLower == jLower && iUpper == jUpper)
    {
      chunk_add_identity(alpha, chunk);
    }
    norm_2 = chunk_get_norm(chunk);
    norm = sqrt(norm_2);
  }

  contribute(cb);
}

#include "node.def.h"