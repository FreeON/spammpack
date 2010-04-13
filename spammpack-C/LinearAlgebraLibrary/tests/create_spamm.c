#include <spamm.h>
#include <stdlib.h>

int
main ()
{
  int M = 10;
  int N = 10;
  int i, j;
  double *A_dense;
  struct spamm_t A;

  A_dense = (double*) malloc(sizeof(double)*M*N);
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      A_dense[spamm_dense_index(i, j, N)] = rand();
    }
  }

  spamm_new(&A);
  spamm_dense_to_spamm(A_dense, &A);
  spamm_print_dense(A_dense, M, N);
  spamm_print_spamm(&A);
}
