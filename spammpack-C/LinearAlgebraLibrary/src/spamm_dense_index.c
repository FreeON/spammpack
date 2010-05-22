/** @file */

#include "spamm.h"

/** Computes a 1-D index from a 2-D index pair.
 *
 * All dense matrices are stored in 1-D NxN vectors. spamm_dense_index() maps
 * a 2-D index into a 1-D index.
 */
int
spamm_dense_index (const int i, const int j, const int M, const int N)
{
  /* Row-Major format. */
  //return i*N+j;

  /* Column-Major format.
   *
   * When using dgemm, we need column-major storage.
   */
  return i+j*M;
}
