#include "spamm.h"

/** Return a linear offset into a dense matrix block in row major order.
 *
 * @param i The row index.
 * @param j The column index.
 * @param M The number of rows of the matrix.
 * @param N The number of columns of the matrix.
 *
 * @return The linear offset into the matrix.
 */
unsigned int
spamm_index_row_major (const unsigned int i, const unsigned int j,
    const unsigned int M, const unsigned int N)
{
  return i*N+j;
}

/** Return a linear offset into a dense matrix block at the kernel tier.
 *
 * @param i The row index.
 * @param j The column index.
 *
 * @return The offset into a kernel tier matrix block.
 */
unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j)
{
  unsigned int offset;

  offset = SPAMM_N_BLOCK*SPAMM_N_BLOCK*spamm_index_row_major(i/SPAMM_N_BLOCK,
      j/SPAMM_N_BLOCK, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)
    + spamm_index_row_major(i%SPAMM_N_BLOCK, j%SPAMM_N_BLOCK, SPAMM_N_BLOCK, SPAMM_N_BLOCK);

  return offset;
}
