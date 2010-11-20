#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int result = 0;
  unsigned int i, j;
  unsigned int N = 64;
  struct spamm_t *A;
  float *A_dense = (float*) malloc(sizeof(float)*N*N);

  for (i = 0; i < N*N; i++)
  {
    A_dense[i] = rand()/(float) RAND_MAX;
  }

  A = spamm_new(N, N);
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      spamm_set(i, j, A_dense[i*N+j], A);
    }
  }

  printf("A_dense:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %1.2f", A_dense[i*N+j]);
    }
    printf("\n");
  }

  /* For debugging, print out the whole tree. */
  spamm_print(A);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (A_dense[i*N+j] != spamm_get(i, j, A))
      {
        result = -1;
        printf("failed test at A[%u][%u] (found %f, should be %f)\n", i, j,
            spamm_get(i, j, A), A_dense[i*N+j]);
        break;
      }
    }
    if (result < 0) { break; }
  }

  spamm_delete(&A);
  free(A_dense);
  return result;
}
