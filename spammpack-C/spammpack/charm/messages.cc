/** @file
 *
 * The implementation of the various messages.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "messages.h"
#include "logger.h"

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
 * @param nodes The Node array.
 */
MatrixNodeMsg::MatrixNodeMsg (CProxy_Node nodes)
{
  this->nodes = nodes;
}

/** The constructor.
 *
 * @param index The linear index of this node.
 * @param norm The norm of this node.
 * @param norm_2 The square of the norm of this node.
 */
NodeInfoMsg::NodeInfoMsg (int index, double norm, double norm_2)
{
  this->index = index;
  this->norm = norm;
  this->norm_2 = norm_2;
}

#include "messages.def.h"
