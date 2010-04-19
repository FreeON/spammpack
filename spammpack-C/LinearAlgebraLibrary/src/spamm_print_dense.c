#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_print_dense (const int M, const int N, const double *A_dense)
{
  assert(A_dense != NULL);
  assert(M > 0 && N > 0);

  int i, j;

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      printf(" % f", A_dense[spamm_dense_index(i, j, M, N)]);
    }
    printf("\n");
  }
}
