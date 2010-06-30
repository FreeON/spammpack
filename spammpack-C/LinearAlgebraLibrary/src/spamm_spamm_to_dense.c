#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** Convert a spamm_t matrix into a dense matrix.
 *
 * This function allocates the memory for A_dense. The caller needs to take
 * care of free()'ing that memory.
 *
 * @param A The spamm_t matrix.
 * @param A_dense The dense matrix.
 */
void
spamm_spamm_to_dense (const struct spamm_t *A, floating_point_t **A_dense)
{
  int i, j;

  assert(A != NULL);
  assert(A_dense != NULL);

  *A_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*A->M*A->N);

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      (*A_dense)[spamm_dense_index(i, j, A->M, A->N)] = spamm_get(i, j, A);
    }
  }
}
