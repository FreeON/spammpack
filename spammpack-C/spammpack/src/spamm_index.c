#include "spamm.h"
#include <assert.h>

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
  assert(i < M);
  assert(j < N);

  return i*N+j;
}

/** Return a linear offset into a dense matrix block at the kernel tier. The
 * indices are within the kernel block at the kernel tier, i.e. within the
 * range of [0, SPAMM_N_KERNEL[.
 *
 * @param i The row index within the kernel block.
 * @param j The column index within the kernel block.
 *
 * @return The offset into a kernel tier matrix block.
 */
unsigned int
spamm_index_kernel_block (const unsigned int i, const unsigned int j)
{
  unsigned int offset;

  assert(i < SPAMM_N_KERNEL);
  assert(j < SPAMM_N_KERNEL);

  offset = SPAMM_N_BLOCK*SPAMM_N_BLOCK
    * spamm_index_row_major(i/SPAMM_N_BLOCK, j/SPAMM_N_BLOCK, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)
    + spamm_index_row_major(i%SPAMM_N_BLOCK, j%SPAMM_N_BLOCK, SPAMM_N_BLOCK, SPAMM_N_BLOCK);

  return offset;
}

/** Return a linear offset into the norms at the kernel tier. The indices are
 * within the kernel block at the kernel tier, i.e. within the range of [0,
 * SPAMM_N_KERNEL[.
 *
 * @param i The row index within the kernel block.
 * @param j The column index within the kernel block.
 *
 * @return The offset into the norm[] arry at the kernel tier.
 */
unsigned int
spamm_index_norm (const unsigned int i, const unsigned int j)
{
  unsigned int offset;

  assert(i < SPAMM_N_KERNEL);
  assert(j < SPAMM_N_KERNEL);

  offset = spamm_index_row_major(i/SPAMM_N_BLOCK, j/SPAMM_N_BLOCK, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK);

  return offset;
}
