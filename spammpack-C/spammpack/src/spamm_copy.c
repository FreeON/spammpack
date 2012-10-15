/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>

/** Interface to spamm_hashed_delete(). */
void
spamm_hashed_delete (struct spamm_hashed_t **A);

/** Copy a spamm_hashed_data_t object. \f$ data_A \leftarrow \beta data_B \f$.
 *
 * @param beta The scalar beta.
 * @param data_B The matrix to copy from.
 *
 * @return The newly created spamm_data_t object.
 */
struct spamm_hashed_data_t *
spamm_hashed_data_copy (const float beta,
    const struct spamm_hashed_data_t *const data_B)
{
  unsigned int i;
  struct spamm_hashed_data_t *data;

  data = spamm_hashed_new_data(data_B->tier, data_B->index_2D, data_B->layout);

  data->node_norm = beta*data_B->node_norm;
  data->node_norm2 = data->node_norm*data->node_norm;

  for (i = 0; i < SPAMM_N_KERNEL_BLOCKED*SPAMM_N_KERNEL_BLOCKED; i++)
  {
    data->norm[i] = beta*data_B->norm[i];
    data->norm2[i] = data->norm[i]*data->norm[i];
  }

  for (i = 0; i < SPAMM_N_KERNEL*SPAMM_N_KERNEL; i++)
  {
    data->block_dense[i] = beta*data_B->block_dense[i];
    data->block_dense_store[i] = beta*data_B->block_dense_store[i];
    data->block_dense_transpose[i] = beta*data_B->block_dense_transpose[i];
    data->block_dense_dilated[4*i+0] = beta*data_B->block_dense_dilated[4*i+0];
    data->block_dense_dilated[4*i+1] = beta*data_B->block_dense_dilated[4*i+1];
    data->block_dense_dilated[4*i+2] = beta*data_B->block_dense_dilated[4*i+2];
    data->block_dense_dilated[4*i+3] = beta*data_B->block_dense_dilated[4*i+3];
  }

  return data;
}

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
  unsigned int i;
  struct spamm_hashtable_t *tier_hashtable;
  struct spamm_list_t *index;
  struct spamm_hashed_data_t *data;

  assert(B != NULL);

  if (*A != NULL)
  {
    /* Delete what's there already. */
    spamm_hashed_delete(A);
  }

  *A = spamm_hashed_new(B->tier, B->kernel_tier, B->depth,
      B->N_padded,
      B->M_lower, B->M_upper,
      B->N_lower, B->N_upper);

  tier_hashtable = B->tier_hashtable[B->kernel_tier-B->tier];
  index = spamm_hashtable_keys(tier_hashtable);

  if (spamm_list_length(index) != spamm_hashtable_get_number_keys(tier_hashtable))
  {
    SPAMM_FATAL("inconsistent number of keys (list = %u, hashtable = %u)\n",
        spamm_list_length(index), spamm_hashtable_get_number_keys(tier_hashtable));
  }

  for (i = 0; i < spamm_hashtable_get_number_keys(tier_hashtable); i++)
  {
    data = spamm_hashtable_lookup(tier_hashtable, spamm_list_get_index(index, i));
    spamm_hashtable_insert((*A)->tier_hashtable[(*A)->kernel_tier-(*A)->tier], spamm_list_get_index(index, i), spamm_hashed_data_copy(beta, data));
  }
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

  if (B == NULL && *A != NULL)
  {
    spamm_recursive_delete(A);
  }

  if (B != NULL)
  {
    if (*A == NULL)
    {
      *A = spamm_recursive_new_node(B->tier,
          B->number_dimensions,
          B->N_contiguous, B->N_linear,
          B->N_lower, B->N_upper);
    }

    if ((*A)->N_upper[0]-(*A)->N_lower[0] == (*A)->N_linear)
    {
      spamm_hashed_copy(&(*A)->hashed_tree, beta, B->hashed_tree);
    }

    else if ((*A)->N_upper[0]-(*A)->N_lower[0] == (*A)->N_contiguous)
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
  unsigned int i;

  assert(B != NULL);

  if (*A == NULL)
  {
    /* Create new matrix A. */
    *A = spamm_new(B->number_dimensions, B->N, B->linear_tier, B->contiguous_tier, B->layout);
  }

  /* Sanity check. */
  if ((*A)->number_dimensions != B->number_dimensions)
  {
    SPAMM_FATAL("mismatch in number of dimensions\n");
  }

  for (i = 0; i < B->number_dimensions; i++)
  {
    if ((*A)->N[i] != B->N[i])
    {
      SPAMM_FATAL("mismatch in matrix dimension N[%u]\n", i);
    }
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
