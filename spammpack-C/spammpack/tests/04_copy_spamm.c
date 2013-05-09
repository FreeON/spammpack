#include "test.h"

#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

#define REL_TOLERANCE 1e-6
#define ABS_TOLERANCE 1e-8
#define DECAY 2.5

int
main (int argc, char **argv)
{
  int result = 0;

  float beta = 1.2;

  unsigned int number_dimensions;

  unsigned int i;
  unsigned int *N;

  unsigned int N_contiguous;

  unsigned int dim;

  const unsigned int chunk_tier = 5;

  short use_linear_tree;
  enum matrix_t matrix_type;

  struct spamm_matrix_t *A = NULL;
  struct spamm_matrix_t *B = NULL;

  double flop;
  double memop;

  float *B_dense;

  for(number_dimensions = 1; number_dimensions <= 3; number_dimensions++)
  {
    for(use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
    {
      if(number_dimensions != 2 && use_linear_tree)
      {
        continue;
      }

      for(matrix_type = 0; matrix_type < NUMBER_MATRIX_TYPES; matrix_type++)
      {
        printf("dim: %u, linTree: %u, matrix_type: %s, ", number_dimensions, use_linear_tree, get_matrix_type_name(matrix_type));

        if(matrix_type == exponential_decay)
        {
          N = generate_shape(number_dimensions, 1);
        }

        else
        {
          N = generate_shape(number_dimensions, 0);
        }

        N_contiguous = 1;
        for(dim = 0; dim < number_dimensions; dim++)
        {
          N_contiguous *= N[dim];
        }
        B_dense = generate_matrix_float(number_dimensions, matrix_type, N, DECAY);
        B = spamm_convert_dense_to_spamm(number_dimensions, N, chunk_tier, use_linear_tree, row_major, B_dense);

        printf("B info: ");
        spamm_print_info(B);
        if(spamm_check(B, REL_TOLERANCE) != SPAMM_OK)
        {
          SPAMM_FATAL("failed\n");
        }

        spamm_copy(&A, beta, B, &flop, &memop);

        if(spamm_check(A, REL_TOLERANCE) != SPAMM_OK)
        {
          SPAMM_FATAL("failed\n");
        }

        for(i = 0; i < N_contiguous; i++)
        {
          B_dense[i] *= beta;
        }

        if(compare_spamm_to_dense_float(A, B_dense, ABS_TOLERANCE) != SPAMM_OK)
        {
          SPAMM_FATAL("comparison failed\n");
        }

        free(N);
        free(B_dense);
        spamm_delete(&A);
        spamm_delete(&B);
      }
    }
  }

  return result;
}
