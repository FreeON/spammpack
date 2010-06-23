#include "spamm.h"

/** Convert a linear index to a coordinate pair, i.e. index --> (i, j).
 *
 * @param index The linear index.
 * @param i The row index.
 * @param j The column index.
 * @param M The number of rows of the matrix (the padded size).
 * @param N The number of columns of the matrix (the padded size).
 * @param M_block The number of rows of the matrix blocks.
 * @param N_block The number of columns of the matrix blocks.
 */
void
spamm_linear_to_coordinates (const unsigned int index, unsigned int *i,
    unsigned int *j, const unsigned int M, const unsigned int N,
    const unsigned int M_block, const unsigned int N_block)
{
  unsigned int bitmask = 1;
  unsigned int M_lower, N_lower;
  unsigned int M_upper, N_upper;

  M_lower = 0;
  N_lower = 0;
  M_upper = M;
  N_upper = N;
}
