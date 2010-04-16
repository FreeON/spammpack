#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

int
main ()
{
  struct spamm_t A, B, B_test;

  int max_size = 5;
  int M[5] = { 2, 10, 15, 400, 600 };
  int N[5] = { 2, 10,  8, 400, 300 };

  int max_block = 11;
  int M_block[11] = { 1, 1, 1, 2, 2, 2, 5, 8, 8, 10, 100 };
  int N_block[11] = { 1, 2, 3, 1, 2, 5, 2, 8, 2, 10, 100 };

  int max_child = 5;
  int M_child[5] = { 2, 2, 2, 3, 4 };
  int N_child[5] = { 2, 3, 4, 3, 4 };

  int max_fill = 4;
  double fill[4] = { 0.1, 0.2, 0.5, 1.0 };

  int i, j;
  int i_size;
  int i_block;
  int i_child;
  int i_fill;

  for (i_fill = 0; i_fill < max_fill; ++i_fill) {
    for (i_size = 0; i_size < max_size; ++i_size) {
      for (i_block = 0; i_block < max_block; ++i_block) {
        for (i_child = 0; i_child < max_child; ++i_child)
        {
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &A);
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &B);
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &B_test);

          printf("[add_spamm] %ix%i matrix, %.1f%% full, padded %ix%i, i_block dimensions %ix%i, i_child dimensions %ix%i\n",
              M[i_size], N[i_size], fill[i_fill]*100, A.M_padded, A.N_padded, M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child]);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill[i_fill]) ? rand()/(double) RAND_MAX : 0.0), &A);
              spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill[i_fill]) ? rand()/(double) RAND_MAX : 0.0), &B);
            }
          }

          //printf("A =\n");
          //spamm_print_spamm(&A);
          //printf("B =\n");
          //spamm_print_spamm(&B);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              spamm_set(i, j, spamm_get(i, j, &A)+spamm_get(i, j, &B), &B_test);
            }
          }

          //printf("B_test =\n");
          //spamm_print_spamm(&B_test);

          spamm_add(1.0, &A, 1.0, &B);

          //printf("B =\n");
          //spamm_print_spamm(&B);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              if (spamm_get(i, j, &B) != spamm_get(i, j, &B_test))
              {
                printf("[add_spamm] mismatch: (B_test[%i][%i] = %e) != (B[%i][%i] = %e)\n", i, j, spamm_get(i, j, &B_test), i, j, spamm_get(i, j, &B));
                exit(1);
              }
            }
          }

          spamm_delete(&A);
          spamm_delete(&B);
          spamm_delete(&B_test);
        }
      }
    }
  }

  return 0;
}
