#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Convert a dense matrix into a SpAMM matrix.
 *
 * If the spamm_t matrix A is already initialized, then it is checked whether
 * its parameters are identical to the parameters requested. If not, the
 * matrix is deleted and re-initialized.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param A_dense The dense input matrix.
 * @param A The spamm_t output matrix.
 */
void
spamm_dense_to_spamm (const unsigned int M, const unsigned int N,
    const floating_point_t *A_dense, struct spamm_t *A)
{
  int i, j;

  assert(A_dense != NULL);
  assert(A != NULL);

  if (M != A->M || N != A->N)
  {
    /* De- and re-allocate A with the correct dimensions. */
    spamm_delete(A);
    spamm_new(M, N, A);
  }

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      spamm_set(i, j, A_dense[spamm_dense_index(i, j, M, N)], A);
    }
  }
}
