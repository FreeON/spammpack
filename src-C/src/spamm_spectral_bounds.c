/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <math.h>

/** Calculate the spectral bounds on A using the Gershgorin circle theorem.
 *
 * @param [out] a_min The lower spectral bound.
 * @param [out] a_max The upper spectral bound.
 * @param A The matrix.
 */
void
spamm_spectral_bounds (float *const a_min,
    float *const a_max,
    struct spamm_matrix_t *A)
{
  unsigned int i[2];
  float rowsum;
  float a_min_temp;
  float a_max_temp;

  if(A->number_dimensions != 2)
  {
    SPAMM_FATAL("Spectral bounds can only be calculated using this method for rank 2 matrix\n");
  }

  if(A->N[0] != A->N[1])
  {
    SPAMM_FATAL("Spectral bounds can only be calculated for a square matrix\n");
  }

  for(i[0] = 0; i[0] < A->N[0]; i[0]++) {
    for(i[1] = 0, rowsum = 0; i[1] < A->N[1]; i[1]++)
    {
      if(i[0] != i[1])
      {
        rowsum += fabsf(spamm_get(i, A));
      }
    }

    i[1] = i[0];
    a_min_temp = spamm_get(i, A)-rowsum;
    a_max_temp = spamm_get(i, A)+rowsum;

    if(i[0] == 0)
    {
      *a_min = a_min_temp;
      *a_max = a_max_temp;
    }

    else
    {
      if(a_min_temp < *a_min) { *a_min = a_min_temp; }
      if(a_max_temp > *a_max) { *a_max = a_max_temp; }
    }
  }
}
