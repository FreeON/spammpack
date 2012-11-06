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

  if (number_dimensions == 2 && use_linear_tree && tier == chunk_tier)
  {
    spamm_hashed_copy(&(*A)->tree.hashed_tree, beta, B->tree.hashed_tree);
  }

  else if (tier == chunk_tier)
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

  if (B->number_dimensions == 2 && B->use_linear_tree && B->chunk_tier == 0)
  {
    spamm_hashed_copy(&(*A)->tree.hashed_tree, beta, B->tree.hashed_tree);
  }

  else if (B->chunk_tier == 0)
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
