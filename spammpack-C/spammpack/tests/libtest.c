/** @file
 *
 * Generate matrices for the tests.
 */

#include "test.h"

#include <math.h>
#include <spamm.h>
#include <stdlib.h>

/** Generate a matrix shape.
 *
 * @param number_dimensions The number of dimensions.
 * @param is_square Whether the matrix is square shaped.
 *
 * @return The shape of the matrix.
 */
unsigned int *
generate_shape (const unsigned int number_dimensions,
    const short is_square)
{
  unsigned int *N;
  unsigned int dim;

  N = calloc(number_dimensions, sizeof(unsigned int));

  for(dim = 0; dim < number_dimensions; dim++)
  {
    N[dim] = 150+(int) ((0.5-(float) rand()/(float) RAND_MAX)*30);
  }

  return N;
}

/** Generate a random test matrix.
 *
 * @param number_dimensions The number of dimensions of the matrix.
 * @param is_sparse Whether the matrix is sparse or not.
 * @param N The shape of the matrix.
 *
 * @return The newly allocated matrix.
 */
float *
generate_matrix (const unsigned int number_dimensions,
    const short is_sparse,
    const unsigned int *const N)
{
  float *A;

  unsigned int dim;
  unsigned int *i;
  unsigned int N_contiguous;

  i = calloc(number_dimensions, sizeof(unsigned int));

  N_contiguous = 1;
  for(dim = 0; dim < number_dimensions; dim++)
  {
    N_contiguous *= N[dim];
  }

  A = (float*) calloc(N_contiguous, sizeof(float));

  if(is_sparse)
  {
    switch(number_dimensions)
    {
      case 1:
        for(i[0] = 0; i[0] < N[0]; i[0]++)
        {
          A[i[0]] = rand()/(float) RAND_MAX;
        }
        break;

      case 2:
        for(i[0] = 0; i[0] < N[0]; i[0]++) {
          for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++)
          {
            A[spamm_index_row_major(i[0], i[1], N[0], N[1])] = rand()/(float) RAND_MAX;
          }
        }
        break;

      case 3:
        for(i[0] = 0; i[0] < N[0]; i[0]++) {
          for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++) {
            for((i[0] >= 10 ? i[2] = i[0]-10 : 0); i[2] < i[0]+10 && i[2] < N[2]; i[2]++)
            {
              A[i[0]+N[0]*(i[1]+N[1]*i[2])] = rand()/(float) RAND_MAX;
            }
          }
        }
        break;

      default:
        SPAMM_FATAL("FIXME\n");
        break;
    }
  }

  else
  {
    for(i[0] = 0; i[0] < N_contiguous; i[0]++)
    {
      A[i[0]] = rand()/(float) RAND_MAX;
    }
  }

  free(i);

  return A;
}

/** Compare a SpAMM matrix with a dense matrix.
 *
 * @param A The SpAMM matrix.
 * @param A_dense The dense matrix.
 * @param abs_tolerance The absolute tolerance.
 *
 * @return SPAMM_OK in case the matrices are identical, and SPAMM_ERROR in
 */
int
compare_spamm_to_dense (const struct spamm_matrix_t *const A,
    const float *const A_dense,
    const double abs_tolerance)
{
  int result = SPAMM_OK;
  double max_diff;
  unsigned int *max_diff_i;

  unsigned int *i;
  unsigned int *N;

  i = calloc(spamm_get_number_dimensions(A), sizeof(unsigned int));
  N = spamm_get_N(A);

  max_diff = 0.0;
  max_diff_i = calloc(spamm_get_number_dimensions(A), sizeof(unsigned int));

  switch(spamm_get_number_dimensions(A))
  {
    case 1:
      for(i[0] = 0; i[0] < N[0]; i[0]++)
      {
        if(fabs(A_dense[i[0]]-spamm_get(i, A)) > max_diff)
        {
          max_diff = fabs(A_dense[i[0]]-spamm_get(i, A));
          max_diff_i[0] = i[0];
        }
      }

      if(max_diff > abs_tolerance)
      {
        result |= SPAMM_ERROR;
        SPAMM_WARN("max diff of A[%u] (absolute tolerance was %1.2e) = %e\n",
            max_diff_i[0], abs_tolerance, max_diff);
      }
      break;

    case 2:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++)
        {
          if(fabs(A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]-spamm_get(i, A)) > max_diff)
          {
            max_diff = fabs(A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])]-spamm_get(i, A));
            max_diff_i[0] = i[0];
            max_diff_i[1] = i[1];
          }
        }
      }

      if(max_diff > abs_tolerance)
      {
        result |= SPAMM_ERROR;
        SPAMM_WARN("max diff of A[%u][%u] (absolute tolerance was %1.2e) = %e\n",
            max_diff_i[0], max_diff_i[1], abs_tolerance, max_diff);
      }
      break;

    case 3:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++) {
          for(i[2] = 0; i[2] < N[2]; i[2]++)
          {
            if(fabs(A_dense[spamm_index_row_major_3(spamm_get_number_dimensions(A), N, i)]-spamm_get(i, A)) > max_diff)
            {
              max_diff = fabs(A_dense[spamm_index_row_major_3(spamm_get_number_dimensions(A), N, i)]-spamm_get(i, A));
              max_diff_i[0] = i[0];
              max_diff_i[1] = i[1];
              max_diff_i[2] = i[2];
            }
          }
        }
      }

      if(max_diff > abs_tolerance)
      {
        result |= SPAMM_ERROR;
        SPAMM_WARN("max diff of A[%u][%u][%u] (ref = %e, SpAMM = %e, absolute tolerance was %1.2e) = %e\n",
            max_diff_i[0], max_diff_i[1], max_diff_i[2],
            A_dense[spamm_index_row_major_3(spamm_get_number_dimensions(A), N, max_diff_i)],
            spamm_get(max_diff_i, A), abs_tolerance, max_diff);
      }
      break;

    default:
      SPAMM_FATAL("FIXME\n");
      break;
  }

  free(i);

  return SPAMM_OK;
}
