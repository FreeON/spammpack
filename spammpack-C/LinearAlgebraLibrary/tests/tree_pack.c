#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int M = 6;
  int N = 6;

  int M_block = 2;
  int N_block = 1;

  int linear_tier = 0;

  int chunksize = 100;

  int i, j;
  struct spamm_t A;

  spamm_new(M, N, M_block, N_block, 2, 2, 0.0, &A);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      spamm_set(i, j, spamm_dense_index(i, j, M, N), &A);
    }
  }

  printf("A (tree)\n");
  spamm_print_tree(&A);

  printf("A =\n");
  spamm_print_spamm(&A);

  spamm_tree_pack(linear_tier, chunksize, i_mask, &A);

  //printf("A (tree)\n");
  //spamm_print_tree(&A);

  //printf("A =\n");
  //spamm_print_spamm(&A);

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      if (spamm_get(i, j, &A) != spamm_dense_index(i, j, M, N))
      {
        LOG("element mismatch A(%i,%i) = %f but should be %u\n", i, j, spamm_get(i, j, &A), spamm_dense_index(i, j, M, N));
        return -1;
      }
    }
  }

  spamm_delete(&A);

  return 0;
}
