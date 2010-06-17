#include <spamm.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define THRESHOLD 1e-14

int
main ()
{
  struct spamm_t A, B, C, C_test;

  int max_size = 4;
  int M[4] = { 4, 10, 15, 100 };
  int N[4] = { 4, 10, 15, 100 };

  int max_block = 4;
  int M_block[4] = { 1, 2, 10, 25 };
  int N_block[4] = { 1, 2, 10, 25 };

  int max_child = 3;
  int M_child[3] = { 2, 3, 4 };
  int N_child[3] = { 2, 3, 4 };

  int max_fill = 4;
  double fill[4] = { 0.01, 0.2, 0.5, 1.0 };

  double max_diff;
  int max_diff_i = 0;
  int max_diff_j = 0;

  int i, j, k;
  int i_size;
  int i_block;
  int i_child;
  int i_fill;

  int result = 0;

  for (i_fill = 0; i_fill < max_fill; ++i_fill) {
    for (i_size = 0; i_size < max_size; ++i_size) {
      for (i_block = 0; i_block < max_block; ++i_block) {
        for (i_child = 0; i_child < max_child; ++i_child)
        {
          //i_fill = 3;
          //i_size = 0;
          //i_block = 0;
          //i_child = 0;

          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &A);
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &B);
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &C);
          spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &C_test);

          //printf("[multiply_spamm] %ix%i matrix, %.1f%% full, padded %ix%i, i_block dimensions %ix%i, i_child dimensions %ix%i\n",
          //    M[i_size], N[i_size], fill[i_fill]*100, A.M_padded, A.N_padded, M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child]);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill[i_fill]) ? rand()/(double) RAND_MAX : 0.0), &A);
              spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill[i_fill]) ? rand()/(double) RAND_MAX : 0.0), &B);
            }
          }

          //printf("A =\n");
          //spamm_print_spamm(&A);
          //spamm_print_tree(&A);
          //printf("B =\n");
          //spamm_print_spamm(&B);
          //spamm_print_tree(&B);

          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j) {
              for (k = 0; k < M[i_size]; ++k)
              {
                spamm_set(i, j, spamm_get(i, k, &A)*spamm_get(k, j, &B)+spamm_get(i, j, &C_test), &C_test);
              }
            }
          }

          //printf("C_test =\n");
          //spamm_print_spamm(&C_test);

          spamm_multiply(tree, 1.0, &A, &B, 1.0, &C);

          //printf("C =\n");
          //spamm_print_spamm(&C);

          max_diff = 0;
          for (i = 0; i < M[i_size]; ++i) {
            for (j = 0; j < N[i_size]; ++j)
            {
              if (fabs(spamm_get(i, j, &C)-spamm_get(i, j, &C_test)) > max_diff)
              {
                max_diff = fabs(spamm_get(i, j, &C)-spamm_get(i, j, &C_test));
                max_diff_i = i;
                max_diff_j = j;
              }
            }
          }

          if (max_diff > THRESHOLD)
          {
            printf("[multiply_spamm] biggest mismatch above threshold of %e: (C_test[%i][%i] = %e) != (C[%i][%i] = %e), |diff| = %e\n",
                THRESHOLD,
                max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &C_test),
                max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &C),
                fabs(spamm_get(max_diff_i, max_diff_j, &C_test)-spamm_get(max_diff_i, max_diff_j, &C)));
            result = 1;
          }

          spamm_delete(&A);
          spamm_delete(&B);
          spamm_delete(&C);
          spamm_delete(&C_test);
        }
      }
    }
  }

  return result;
}
