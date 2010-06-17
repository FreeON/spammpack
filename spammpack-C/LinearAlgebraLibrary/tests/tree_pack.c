#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  int number_sizes = 4;
  int M[] = { 2, 6, 10, 100 };
  int N[] = { 2, 6, 6,   85 };

  int number_blocks = 2;
  int M_block[] = { 1, 2 };
  int N_block[] = { 1, 1 };

  int number_linear_tiers = 3;
  int linear_tier[] = { 0, 1, 3 };

  int number_chunks = 2;
  int chunksize[] = { 100, 1000 };

  int i, j;
  int i_size, i_block, i_linear, i_chunk;
  struct spamm_t A;

  for (i_size = 0; i_size < number_sizes; ++i_size) {
    for (i_block = 0; i_block < number_blocks; ++i_block) {
      for (i_linear = 0; i_linear < number_linear_tiers; ++i_linear) {
        for (i_chunk = 0; i_chunk < number_chunks; ++i_chunk)
        {
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], 2, 2, 0.0, &A);
          LOG("%ux%u matrix, %ux%u blocks, linear tier = %u, chunksize = %u bytes, depth = %u\n",
              M[i_size], N[i_size], M_block[i_block], N_block[i_block], linear_tier[i_linear], chunksize[i_chunk], A.tree_depth);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              spamm_set(i, j, spamm_dense_index(i, j, M[i_size], N[i_size]), &A);
            }
          }

          //printf("A (tree)\n");
          //spamm_print_tree(&A);

          //printf("A =\n");
          //spamm_print_spamm(&A);

          spamm_tree_pack(linear_tier[i_linear], chunksize[i_chunk], i_mask, &A);

          //printf("A (tree)\n");
          //spamm_print_tree(&A);

          //printf("A =\n");
          //spamm_print_spamm(&A);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              if (spamm_get(i, j, &A) != spamm_dense_index(i, j, M[i_size], N[i_size]))
              {
                LOG("element mismatch A(%i,%i) = %f but should be %u\n", i, j, spamm_get(i, j, &A), spamm_dense_index(i, j, M[i_size], N[i_size]));
                return -1;
              }
            }
          }

          spamm_delete(&A);
          //exit(0);
        }
      }
    }
  }

  return 0;
}
