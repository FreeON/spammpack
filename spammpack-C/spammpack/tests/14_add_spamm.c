#include "test.h"

#include <spamm.h>
#include <math.h>
#include <stdlib.h>

#define REL_TOLERANCE 1e-8
#define TEST_TOLERANCE 5e-7
#define DECAY 2.5

int
main (int argc, char **argv)
{
  int result = 0;

  unsigned int number_dimensions;

  unsigned int *N;
  unsigned int *i;

  const unsigned int chunk_tier = 4;

  short use_linear_tree;
  short symmetric;

  enum matrix_t matrix_type_A;
  enum matrix_t matrix_type_B;

  double alpha = 1.2;
  double beta = 0.8;

  float *A_dense;
  float *B_dense;

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;

  double flop;
  double memop;

  for(number_dimensions = 1; number_dimensions <= 3; number_dimensions++) {
    for(use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
    {
      if(number_dimensions != 2 && use_linear_tree)
      {
        continue;
      }

      for(matrix_type_A = 0; matrix_type_A < NUMBER_MATRIX_TYPES; matrix_type_A++) {
        for(matrix_type_B = 0; matrix_type_B < NUMBER_MATRIX_TYPES; matrix_type_B++)
        {
          SPAMM_INFO("dim: %u, linTree: %u, matrix_type_A: %s, matrix_type_B: %s\n",
              number_dimensions, use_linear_tree,
              get_matrix_type_name(matrix_type_A), get_matrix_type_name(matrix_type_B));

          if(matrix_type_A == exponential_decay || matrix_type_B == exponential_decay)
          {
            symmetric = 1;
          }

          else
          {
            symmetric = 0;
          }

          i = calloc(number_dimensions, sizeof(unsigned int));
          N = generate_shape(number_dimensions, symmetric);

          A_dense = generate_matrix_float(number_dimensions, matrix_type_A, N, DECAY);
          B_dense = generate_matrix_float(number_dimensions, matrix_type_B, N, DECAY);

          /* Convert to SpAMM matrix. */
          A = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, A_dense);
          B = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, B_dense);

          SPAMM_INFO("A: ");
          spamm_print_info(A);

          SPAMM_INFO("B: ");
          spamm_print_info(B);

          /* Add by hand for verification. */
          switch(number_dimensions)
          {
            case 1:
              for(i[0] = 0; i[0] < N[0]; i[0]++)
              {
                A_dense[i[0]] = alpha*A_dense[i[0]] + beta*B_dense[i[0]];
              }
              break;

            case 2:
              for(i[0] = 0; i[0] < N[0]; i[0]++) {
                for(i[1] = 0; i[1] < N[1]; i[1]++)
                {
                  A_dense[spamm_index_row_major_3(number_dimensions, N, i)] =
                    alpha*A_dense[spamm_index_row_major_3(number_dimensions, N, i)]
                    + beta*B_dense[spamm_index_row_major_3(number_dimensions, N, i)];
                }
              }
              break;

            case 3:
              for(i[0] = 0; i[0] < N[0]; i[0]++) {
                for(i[1] = 0; i[1] < N[1]; i[1]++) {
                  for(i[2] = 0; i[2] < N[2]; i[2]++)
                  {
                    A_dense[spamm_index_row_major_3(number_dimensions, N, i)] =
                      alpha*A_dense[spamm_index_row_major_3(number_dimensions, N, i)]
                      + beta*B_dense[spamm_index_row_major_3(number_dimensions, N, i)];
                  }
                }
              }
              break;

            default:
              SPAMM_FATAL("FIXME\n");
              break;
          }

          flop = 0;
          memop = 0;
          spamm_add(alpha, A, beta, B, &flop, &memop);
          SPAMM_INFO("%e flop, %e memop\n", flop, memop);

          /* Check tree consistency. */
          if(spamm_check(A, REL_TOLERANCE) != SPAMM_OK)
          {
            SPAMM_FATAL("failed\n");
          }

          if(compare_spamm_to_dense_float(A, A_dense, TEST_TOLERANCE) != SPAMM_OK)
          {
            SPAMM_FATAL("comparison failed\n");
          }

          free(i);
          free(N);
          free(A_dense);
          free(B_dense);
          spamm_delete(&A);
          spamm_delete(&B);
        }
      }
    }
  }

  return result;
}
