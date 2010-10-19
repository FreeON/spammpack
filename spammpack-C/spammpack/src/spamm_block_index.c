#include "spamm.h"

/** Computes a 1-D index from a 2-D index pair.
 *
 * This function calls spamm_dense_index() on the individual blocks to be able
 * to address matrix elements in a contiguous dense block at or below the
 * kernel_tier.
 *
 * @param i The row index.
 * @param j The column index.
 * @param M_block The number of rows of the matrix blocks.
 * @param N_block The number of columns of the matrix blocks.
 * @param M_kernel The number of rows of basic matrix blocks.
 * @param N_kernel The number of columns of basic matrix blocks.
 */
unsigned int
spamm_block_index (const unsigned int i, const unsigned int j,
    const unsigned int M_block, const unsigned int N_block,
    const unsigned int M_kernel, const unsigned int N_kernel)
{
  /* Store blocks... */
}
