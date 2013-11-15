/** @file
 *
 * The implementation of the chunk functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "chunk.h"

struct chunk_t
{
  /** The size of this matrix block. */
  int blocksize;

  /** The matrix data. */
  double *block;
};

/** Get the chunksize.
 */
size_t
chunk_sizeof (const int blocksize)
{
  return sizeof(struct chunk_t)+sizeof(double)*blocksize*blocksize;
}
