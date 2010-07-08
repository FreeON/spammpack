#include "spamm.h"

/** Swap to arrays of type floating_point_t.
 *
 * @param M The number of rows of the dense block.
 * @param N The number of columns of the dense block.
 * @param A The first dense block.
 * @param B The second dense block.
 */
void
spamm_swap_block_dense (const unsigned int M, const unsigned int N, floating_point_t *A, floating_point_t *B)
{
  unsigned int i;
  floating_point_t temp;

  for (i = 0; i < M*N; ++i)
  {
    temp = A[i];
    A[i] = B[i];
    B[i] = temp;
  }
}
