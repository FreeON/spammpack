#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Print out a dense matrix.
 *
 * @param M Number of rows of dense matrix.
 * @param N Number of columns of dense matrix.
 * @param A_dense The dense matrix.
 */
void
spamm_print_dense (const unsigned int M, const unsigned int N, const floating_point_t *A_dense)
{
  int i, j;

  assert(A_dense != NULL);
  assert(M > 0 && N > 0);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      printf(" % f", A_dense[spamm_dense_index(i, j, M, N)]);
    }
    printf("\n");
  }
}
