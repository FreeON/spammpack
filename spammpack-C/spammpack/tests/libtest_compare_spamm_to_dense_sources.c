/** @file */

#include "spamm.h"
#include "test.h"

#include <math.h>

/** Compare a SpAMM matrix with a dense matrix.
 *
 * @param A The SpAMM matrix.
 * @param A_dense The dense matrix.
 * @param abs_tolerance The absolute tolerance.
 *
 * @return SPAMM_OK in case the matrices are identical, and SPAMM_ERROR in
 */
int
SPAMM_FUNCTION(compare_spamm_to_dense, MATRIX_TYPE) (const struct spamm_matrix_t *const A,
    const MATRIX_TYPE *const A_dense,
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
