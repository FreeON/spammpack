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
spamm_recursive_copy (struct spamm_recursive_node_t **A,
    const struct spamm_recursive_node_t *const B)
{
  SPAMM_FATAL("[FIXME]\n");
}

/** Copy a matrix. \f$ A \leftarrow B \f$.
 *
 * @param A The matrix to copy to.
 * @param B The matrix to copy from.
 */
void
spamm_hashed_copy (struct spamm_hashed_t **A,
    const struct spamm_hashed_t *const B)
{
  SPAMM_FATAL("[FIXME]\n");
}

/** Copy a matrix. \f$ A \leftarrow B \f$.
 *
 * @param A The matrix to copy to.
 * @param B The matrix to copy from.
 */
void
spamm_copy (struct spamm_matrix_t **A,
    const struct spamm_matrix_t *const B)
{
  assert(B != NULL);

  if (*A == NULL)
  {
    /* Create new matrix A. */
    *A = spamm_new(B->M, B->N, B->linear_tier, B->contiguous_tier, B->layout);
  }

  /* Sanity check. */
  if ((*A)->M != B->M)
  {
    SPAMM_FATAL("mismatch in matrix dimension M\n");
  }

  if ((*A)->N != B->N)
  {
    SPAMM_FATAL("mismatch in matrix dimension N\n");
  }

  if ((*A)->linear_tier != B->linear_tier)
  {
    SPAMM_FATAL("mismatch in linear tier\n");
  }

  if ((*A)->contiguous_tier != B->contiguous_tier)
  {
    SPAMM_FATAL("mismatch in contiguous tier\n");
  }

  if ((*A)->layout != B->layout)
  {
    SPAMM_FATAL("mismatch in layout\n");
  }

  if (B->recursive_tree != NULL)
  {
    spamm_recursive_copy(&(*A)->recursive_tree, B->recursive_tree);
  }

  else if (B->hashed_tree != NULL)
  {
    spamm_hashed_copy(&(*A)->hashed_tree, B->hashed_tree);
  }
}
