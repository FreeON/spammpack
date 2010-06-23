#include <spamm.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define THRESHOLD 1e-5

int
main ()
{
  struct spamm_t A, B, C, C_test;

  int max_size = 4;
  int M[] = { 4, 10, 15, 30 };
  int N[] = { 4, 10, 15, 20 };

  int max_block = 3;
  int M_block[] = { 1, 2, 10 };
  int N_block[] = { 1, 2, 10 };

  int max_child = 1;
  int M_child[] = { 2 };
  int N_child[] = { 2 };

  int max_fill_A = 4;
  double fill_A[] = { 0.01, 0.2, 0.5, 1.0 };

  int max_fill_B = 4;
  double fill_B[] = { 0.01, 0.2, 0.5, 1.0 };

  int number_linear_tiers = 3;
  int linear_tier[] = { 0, 1, 3 };

  int number_chunks = 3;
  int chunksize[] = { 500, 1000, 10000 };

  double max_diff;
  int max_diff_i = 0;
  int max_diff_j = 0;

  int i, j, k;
  int i_size;
  int i_block;
  int i_child;
  int i_fill_A;
  int i_fill_B;
  int i_linear;
  int i_chunk;

  int result = 0;

  //spamm_set_loglevel(info);

  for (i_fill_A = 0; i_fill_A < max_fill_A; ++i_fill_A) {
    for (i_fill_B = 0; i_fill_B < max_fill_B; ++i_fill_B) {
      for (i_size = 0; i_size < max_size; ++i_size) {
        for (i_block = 0; i_block < max_block; ++i_block) {
          for (i_child = 0; i_child < max_child; ++i_child) {
            for (i_linear = 0; i_linear < number_linear_tiers; ++i_linear) {
              for (i_chunk = 0; i_chunk < number_chunks; ++i_chunk)
              {
                //i_fill_A = 0; i_fill_B = 0; i_size = 3; i_block = 0; i_child = 0; i_linear = 0; i_chunk = 0;

                //printf("i_fill_A = %i; i_fill_B = %i; i_size = %i; i_block = %i; i_child = %i; i_linear = %i; i_chunk = %i;\n",
                //    i_fill_A, i_fill_B, i_size, i_block, i_child, i_linear, i_chunk);

                spamm_new(M[i_size], N[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &A);
                spamm_new(N[i_size], M[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &B);
                spamm_new(M[i_size], M[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &C);
                spamm_new(M[i_size], M[i_size], M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child], 1e-10, &C_test);

                //printf("[multiply_spamm] %ix%i matrix, %.1f%% full, padded %ix%i, i_block dimensions %ix%i, i_child dimensions %ix%i\n",
                //    M[i_size], N[i_size], fill[i_fill]*100, A.M_padded, A.N_padded, M_block[i_block], N_block[i_block], M_child[i_child], N_child[i_child]);

                for (i = 0; i < M[i_size]; ++i) {
                  for (j = 0; j < N[i_size]; ++j)
                  {
                    spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill_A[i_fill_A]) ? rand()/(double) RAND_MAX : 0.0), &A);
                    //spamm_set(i, j, spamm_dense_index(i, j, M[i_size], N[i_size]), &A);
                  }
                }

                for (i = 0; i < N[i_size]; ++i) {
                  for (j = 0; j < M[i_size]; ++j)
                  {
                    spamm_set(i, j, (rand()/(double) RAND_MAX > (1-fill_B[i_fill_B]) ? rand()/(double) RAND_MAX : 0.0), &B);
                    //spamm_set(i, j, spamm_dense_index(i, j, M[i_size], N[i_size]), &B);
                  }
                }

                //printf("A =\n");
                //spamm_print_spamm(&A);
                //spamm_print_tree(&A);
                //printf("B =\n");
                //spamm_print_spamm(&B);
                //spamm_print_tree(&B);

                LOG2_INFO("multiplying to get C_test\n");
                for (i = 0; i < M[i_size]; ++i) {
                  for (j = 0; j < M[i_size]; ++j) {
                    for (k = 0; k < N[i_size]; ++k)
                    {
                      spamm_set(i, j, spamm_get(i, k, &A)*spamm_get(k, j, &B)+spamm_get(i, j, &C_test), &C_test);
                    }
                  }
                }

                //printf("C_test =\n");
                //spamm_print_spamm(&C_test);

                LOG2_INFO("packing trees\n");
                spamm_tree_pack(linear_tier[i_linear], chunksize[i_chunk], i_mask, &A);
                spamm_tree_pack(linear_tier[i_linear], chunksize[i_chunk], i_mask, &B);

                LOG2_INFO("multiplying with spamm\n");
                spamm_multiply(tree, 1.0, &A, &B, 1.0, &C);

                //printf("C =\n");
                //spamm_print_spamm(&C);

                max_diff = 0;
                for (i = 0; i < M[i_size]; ++i) {
                  for (j = 0; j < M[i_size]; ++j)
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
                  printf("i_fill_A = %i; i_fill_B = %i; i_size = %i; i_block = %i; i_child = %i; i_linear = %i; i_chunk = %i;\n",
                      i_fill_A, i_fill_B, i_size, i_block, i_child, i_linear, i_chunk);
                  LOG_FATAL("biggest mismatch above threshold of %e: (C_test[%i][%i] = %e) != (C[%i][%i] = %e), |diff| = %e\n",
                      THRESHOLD,
                      max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &C_test),
                      max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &C),
                      fabs(spamm_get(max_diff_i, max_diff_j, &C_test)-spamm_get(max_diff_i, max_diff_j, &C)));
                  return -1;
                }

                spamm_delete(&A);
                spamm_delete(&B);
                spamm_delete(&C);
                spamm_delete(&C_test);

                //return 0;
              }
            }
          }
        }
      }
    }
  }

  return result;
}
