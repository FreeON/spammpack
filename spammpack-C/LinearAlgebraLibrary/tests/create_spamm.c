#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int M = 10;
  int N = 10;

  int M_block = 1;
  int N_block = 1;

  int i, j;
  double *A_dense;
  struct spamm_t A;

  A_dense = (double*) malloc(sizeof(double)*M*N);
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      A_dense[spamm_dense_index(i, j, N)] = rand()/(double) RAND_MAX;
    }
  }

  spamm_new(M, N, M_block, N_block, &A);
  spamm_dense_to_spamm(M, N, M_block, N_block, A_dense, &A);

  printf("A (dense) =\n");
  spamm_print_dense(M, N, A_dense);

  printf("A (SpAMM) =\n");
  spamm_print_spamm(&A);
}
