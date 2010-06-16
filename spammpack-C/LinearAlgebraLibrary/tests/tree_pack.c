#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int M = 2;
  int N = 2;

  int i, j;
  struct spamm_t A;

  spamm_new(M, N, 1, 1, 2, 2, 0.0, &A);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      //spamm_set(i, j, rand()/(float) RAND_MAX, &A);
      spamm_set(i, j, spamm_dense_index(i, j, M, N), &A);
    }
  }

  printf("A (tree)\n");
  spamm_print_tree(&A);

  printf("A =\n");
  spamm_print_spamm(&A);

  spamm_tree_pack(0, 100, i_mask, &A);

  printf("A (tree)\n");
  spamm_print_tree(&A);

  printf("A =\n");
  spamm_print_spamm(&A);

  return 0;
}
