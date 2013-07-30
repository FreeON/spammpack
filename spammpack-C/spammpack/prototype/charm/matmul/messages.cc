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
 * @param x The double value.
 */
DoubleMsg::DoubleMsg (double x)
{
  this->x = x;
}

/** The constructor.
 *
 * @param i The integer value.
 */
IntMsg::IntMsg (int i)
{
  this->i = i;
}

/** The contructor.
 *
 * @param N The matrix size.
 * @param blocksize The submatrix size at the lowest tier.
 * @param depth The tree depth of the matrix.
 * @param tierNode The array of Node objects. There is one array per tier.
 */
MatrixInfoMsg::MatrixInfoMsg (int N, int blocksize, int depth)
{
  this->N = N;
  this->blocksize = blocksize;
  this->depth = depth;
  this->tierNode = NULL;
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
