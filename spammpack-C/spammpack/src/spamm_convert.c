#include "config.h"
#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Convert a dense matrix to SpAMM.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param linear_tier The linear tier (see spamm_new()).
 * @param contiguous_tier The contiguous tier (see spamm_new()).
 * @param dense_type The storage type of the dense matrix.
 * @param A_dense The dense matrix.
 * @param spamm_layout The layout of the SpAMM data nodes.
 *
 * @return The SpAMM matrix.
 */
struct spamm_matrix_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N,
    const unsigned int linear_tier,
    const unsigned int contiguous_tier,
    const enum spamm_layout_t dense_type,
    const float *const A_dense,
    const enum spamm_layout_t spamm_layout)
{
  struct spamm_matrix_t *A = NULL;
  unsigned int i;
  unsigned int j;

  assert(A_dense != NULL);

  A = spamm_new(M, N, linear_tier, contiguous_tier, spamm_layout);

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++)
    {
      switch(dense_type)
      {
        case row_major:
          spamm_set(i, j, A_dense[spamm_index_row_major(i, j, M, N)], A);
          break;

        case column_major:
          spamm_set(i, j, A_dense[spamm_index_column_major(i, j, M, N)], A);
          break;

        default:
          printf("unknown type\n");
          exit(1);
          break;
      }
    }
  }

  return A;
}
