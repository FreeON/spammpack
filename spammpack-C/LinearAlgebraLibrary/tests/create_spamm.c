#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int M = 1000;
  int N = 1000;

  int M_block[10] = { 1, 1, 1, 2, 2, 2, 5, 8, 8, 10 };
  int N_block[10] = { 1, 2, 3, 1, 2, 5, 2, 8, 2, 10 };

  int i, j, k;
  double *A_dense;
  struct spamm_t A;

  A_dense = (double*) malloc(sizeof(double)*M*N);
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      A_dense[spamm_dense_index(i, j, N)] = rand()/(double) RAND_MAX;
    }
  }

  for (k = 0; k < 10; ++k)
  {
    spamm_new(M, N, M_block[k], N_block[k], &A);
    spamm_dense_to_spamm(M, N, M_block[k], N_block[k], A_dense, &A);

    for (i = 0; i < M; ++i) {
      for (j = 0; j < N; ++j)
      {
        if (A_dense[spamm_dense_index(i, j, N)] != spamm_get(i, j, &A))
        {
          printf("mismatch: A_dense[%i][%i] = %f != A[%i][%i] = %f\n", i, j, A_dense[spamm_dense_index(i, j, N)], i, j, spamm_get(i, j, &A));
          exit(1);
        }
      }
    }

    spamm_delete(&A);
  }

  fprintf(stderr, "test ok\n");
  return 0;
}
