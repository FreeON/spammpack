#include "spamm.h"

/** Computes a 1-D index from a 2-D index pair.
 *
 * All dense matrices are stored in 1-D NxN vectors. spamm_dense_index() maps
 * a 2-D index into a 1-D index.
 *
 * @param i The row index.
 * @param j The column index.
 * @param M The number of rows of the matrix.
 * @param N The number of columns of the matrix.
 *
 * @return The linear index of MxN matrix with index pair i and j.
 */
int
spamm_dense_index (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N)
{
  /* Row-Major format. */
  return i*N+j;
}
