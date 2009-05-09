#include "lal.h"

#include <assert.h>

/* The transpose function returns a newly allocated matrix. It is the caller's
 * responsibility to properly free() this matrix when it is not needed
 * anymore.
 */
lal_matrix_t *
lal_transpose (lal_matrix_t *A)
{
  assert(A->M > 0 && A->N > 0);

  int i, j;
  lal_matrix_t *A_transpose;

  /* Allocate space. */
  lal_allocate(A->N, A->M, &A_transpose);

  for (i = 0; i < A->M; ++i) {
    for (j = 0; j < A->N; ++j)
    {
      lal_set(j, i, lal_get(i, j, A), A_transpose);
    }
  }

  return A_transpose;
}
