/** @file
 *
 * The implementation of the various messages.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "messages.h"

DoubleMsg::DoubleMsg (double x)
{
  this->x = x;
}

IntMsg::IntMsg (int i)
{
  this->i = i;
}

/** The contructor.
 */
MatrixInfoMsg::MatrixInfoMsg (int N, int blocksize, int depth,
    CProxy_Node tierNode)
{
  this->N = N;
  this->blocksize = blocksize;
  this->depth = depth;
  this->tierNode = tierNode;
}

#include "messages.def.h"
