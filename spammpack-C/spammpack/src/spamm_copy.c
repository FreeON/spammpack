/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>

/** Copy a matrix. \f$ A \leftarrow B \f$.
 *
 * @param A The matrix to copy to.
 * @param B The matrix to copy from.
 */
void
spamm_copy (struct spamm_matrix_t *A,
    const struct spamm_matrix_t *const B)
{
  assert(B != NULL);

  if (A == NULL)
  {
    /* Create new matrix A. */
    spamm_new(B->M, B->N, B->linear_tier, B->contiguous_tier, B->layout);
  }
}
