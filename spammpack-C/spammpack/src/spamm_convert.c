#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Convert a dense matrix to SpAMM.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param A_dense The dense matrix.
 *
 * @return The SpAMM matrix.
 */
struct spamm_t *
spamm_convert_dense_to_spamm (const unsigned int M, const unsigned int N, float *A_dense)
{
  struct spamm_t *A = NULL;
  unsigned int i, j;

  assert(A_dense != NULL);

  A = spamm_new(M, N);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      spamm_set(i, j, A_dense[spamm_index_column_major(i, j, N, N)], A);
    }
  }

  return A;
}
