#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

void
spamm_dense_to_spamm (const int M, const int N, const int M_block, const int N_block, const double *A_dense, struct spamm_t *A)
{
  assert(A_dense != NULL);
  assert(A != NULL);

  int i, j;

  if (M != A->M && N != A->N)
  {
    /* De- and re-allocate A with the correct dimensions. */
    spamm_delete(A);
    spamm_new(M, N, M_block, N_block, A);
  }

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      spamm_set(i, j, A_dense[spamm_dense_index(i, j, M)], A);
    }
  }
}
