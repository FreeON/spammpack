/** @file */

#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Convert a dense matrix into a SpAMM matrix.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param M_block Number of rows of matrix block in spamm_t.
 * @param N_block Number of columns of matrix block in spamm_t.
 * @param M_child Number of rows of children array in spamm_node_t.
 * @param N_child Number of columns of children array in spamm_node_t.
 * @param threshold Threshold below which matrix elements are considered zero.
 * @param A_dense The dense input matrix.
 * @param A The spamm_t output matrix.
 */
void
spamm_dense_to_spamm (const int M, const int N, const int M_block,
    const int N_block, const int M_child, const int N_child,
    const float_t threshold, const float_t *A_dense, struct spamm_t *A)
{
  int i, j;

  assert(A_dense != NULL);
  assert(A != NULL);

  if (M != A->M || N != A->N || M_block != A->M_block || N_block != A->N_block || M_child != A->M_child || N_child != A->N_child)
  {
    /* De- and re-allocate A with the correct dimensions. */
    spamm_delete(A);
    spamm_new(M, N, M_block, N_block, M_child, N_child, threshold, A);
  }

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      spamm_set(i, j, A_dense[spamm_dense_index(i, j, M, N)], A);
    }
  }
}
