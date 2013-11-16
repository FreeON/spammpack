/** @file
 *
 * The implementation of the various messages.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "messages.h"
#include "logger.h"

#include <assert.h>

/** The constructor.
 *
 * @param chunksize The size of the chunk.
 * @param chunk The chunk.
 */
ChunkMsg::ChunkMsg (size_t chunksize, char *chunk)
{
  this->chunksize = chunksize;
  if(chunk != NULL)
  {
    DEBUG("copying chunk at %p into chunk at %p\n", chunk, this->chunk);
    memcpy(this->chunk, chunk, chunksize);
  }
}

/** The constructor.
 *
 * @param x The value of this message.
 */
DoubleMsg::DoubleMsg (double x)
{
  this->x = x;
}

/** The constructor.
 *
 * @param i The value of this message.
 */
IntMsg::IntMsg (int i)
{
  this->i = i;
}

/** The constructor.
 *
 * @param N The matrix size.
 * @param blocksize The submatrix size at the lowest tier.
 * @param depth The tree depth of the matrix.
 * @param NPadded The padded matrix size.
 */
MatrixInfoMsg::MatrixInfoMsg (int N, int blocksize, int depth, int NPadded)
{
  this->N = N;
  this->blocksize = blocksize;
  this->depth = depth;
  this->NPadded = NPadded;
}

/** Compare with another MatrixInfoMsg object.
 *
 * @param b The message to compare to.
 *
 * @return Whether they are the same or not.
 */
bool MatrixInfoMsg::equal (MatrixInfoMsg *b)
{
  if(this->N != b->N)
  {
    INFO("mismatch in N\n");
    return false;
  }

  return true;
}

/** The constructor.
 *
 * @param index The linear index of this node.
 * @param iLower The lower row index.
 * @param iUpper The upper row index.
 * @param jLower The lower column index.
 * @param jUpper The upper column index.
 * @param norm The norm of this node.
 * @param norm_2 The square of the norm of this node.
 */
NodeInfoMsg::NodeInfoMsg (int index,
    int iLower,
    int iUpper,
    int jLower,
    int jUpper,
    double norm,
    double norm_2)
{
  this->index = index;
  this->iLower = iLower;
  this->iUpper = iUpper;
  this->jLower = jLower;
  this->jUpper = jUpper;
  this->norm = norm;
  this->norm_2 = norm_2;
}

/** The constructor.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 */
DenseMatrixMsg::DenseMatrixMsg (int M, int N)
{
  memset(A, 0, sizeof(double)*M*N);
}

/** The constructor.
 *
 * @param numberTiers The number of tiers.
 */
MatrixNodeMsg::MatrixNodeMsg (int numberTiers)
{
  this->numberTiers = numberTiers;

  /* We can't use the assignment CProxy::operator=() here since it will try to
   * figure out the delegation status of the lhs, i.e. in nodes[tier] =
   * CProxy_Node::ckNew(), and might segfault if the lhs was not set to 0.
   */
  memset(nodes, 0, numberTiers*sizeof(CProxy_Node));
}

#include "messages.def.h"
