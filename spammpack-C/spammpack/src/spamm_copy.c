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
spamm_recursive_copy (struct spamm_recursive_node_t **A,
    const float beta,
    const struct spamm_recursive_node_t *const B,
    const unsigned int number_dimensions,
    const unsigned int tier,
    const unsigned int chunk_tier,
    const short use_linear_tree)
{
  short i;

  if (*A == NULL)
  {
    *A = spamm_recursive_new_node();
  }

  if (tier == chunk_tier)
  {
    spamm_chunk_copy((*A)->tree.chunk, beta, B->tree.chunk);
  }

  else
  {
    for (i = 0; i < ipow(2, number_dimensions); i++)
    {
      spamm_recursive_copy(&(*A)->tree.child[i], beta, B->tree.child[i],
          number_dimensions, tier+1, chunk_tier, use_linear_tree);
    }
  }
}

/** Copy a matrix. \f$ A \leftarrow \beta B \f$.
 *
 * @param A The matrix to copy to.
 * @param beta The scalar beta.
 * @param B The matrix to copy from.
 */
void
spamm_copy (struct spamm_matrix_t **A,
    const float beta,
    const struct spamm_matrix_t *const B)
{
  spamm_delete(A);

  if (B->chunk_tier == 0)
  {
    spamm_chunk_copy(&(*A)->tree.chunk, beta, B->tree.chunk);
  }

  else
  {
    spamm_recursive_copy(&(*A)->tree.recursive_tree, beta,
        B->tree.recursive_tree, B->number_dimensions, 0, B->chunk_tier,
        B->use_linear_tree);
  }
}
