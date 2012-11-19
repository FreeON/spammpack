#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Check the internal consistency of a matrix.
 *
 * @param A The matrix to check
 * @param tolerance The absolute tolerance when comparing values.
 *
 * @return The following error codes are returned:
 *   - SPAMM_OK - The matrix is consistent.
 *   - SPAMM_ERROR - Something is not consistent.
 */
int
spamm_check (const struct spamm_matrix_t *A, const float tolerance)
{
  int result = 0;

  /* Check norms. */
  //SPAMM_FATAL("FIXME\n");

  return result;
}
