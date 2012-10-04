/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>

/** Copy a matrix. \f$ A \leftarrow \beta B \f$.
 *
 * @param A The matrix to copy to.
 * @param beta The scalar beta.
 * @param B The matrix to copy from.
 */
void
spamm_hashed_copy (struct spamm_hashed_t **A,
    const float beta,
    const struct spamm_hashed_t *const B)
{
  assert(B != NULL);

  if (*A == NULL)
  {
    spamm_hashed_new(B->tier, B->kernel_tier, B->depth,
        B->M_lower, B->M_upper,
        B->N_lower, B->N_upper);
  }

  SPAMM_FATAL("[FIXME]\n");
}

/** Copy a matrix. \f$ A \leftarrow \beta B \f$.
 *
 * @param A The matrix to copy to.
 * @param beta The scalar beta.
 * @param B The matrix to copy from.
 */
void
spamm_recursive_copy (struct spamm_recursive_node_t **A,
    const float beta,
    const struct spamm_recursive_node_t *const B)
{
  unsigned int i;

  assert(B != NULL);

  if (*A == NULL)
  {
    *A = spamm_recursive_new_node(B->tier,
        B->N_contiguous, B->N_linear,
        B->M_lower, B->M_upper,
        B->N_lower, B->N_upper);
  }

  if ((*A)->M_upper-(*A)->M_lower == (*A)->N_linear)
  {
    spamm_hashed_copy(&(*A)->hashed_tree, beta, B->hashed_tree);
  }

  else if ((*A)->M_upper-(*A)->M_lower == (*A)->N_contiguous)
  {
    for (i = 0; i < (*A)->N_contiguous*(*A)->N_contiguous; i++)
    {
      (*A)->data[i] = beta*B->data[i];
      (*A)->norm = beta*B->norm;
      (*A)->norm2 = (*A)->norm*(*A)->norm;
    }
  }

  else
  {
    spamm_recursive_copy(&(*A)->child[0], beta, B->child[0]);
    spamm_recursive_copy(&(*A)->child[1], beta, B->child[1]);
    spamm_recursive_copy(&(*A)->child[2], beta, B->child[2]);
    spamm_recursive_copy(&(*A)->child[3], beta, B->child[3]);
  }
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
    spamm_recursive_copy(&(*A)->recursive_tree, 1.0, B->recursive_tree);
  }

  else if (B->hashed_tree != NULL)
  {
    spamm_hashed_copy(&(*A)->hashed_tree, 1.0, B->hashed_tree);
  }
}
