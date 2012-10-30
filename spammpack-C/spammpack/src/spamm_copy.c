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
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const short use_linear_tree)
{
  unsigned int i;

  float *A_matrix;
  float *B_matrix;

  if (B == NULL && *A != NULL)
  {
    spamm_recursive_delete(number_dimensions, tier, contiguous_tier, use_linear_tree, A);
  }

  if (B != NULL)
  {
    if (*A == NULL)
    {
      *A = spamm_recursive_new_node();
    }

    if (tier == contiguous_tier && use_linear_tree)
    {
      spamm_hashed_copy(&(*A)->tree.hashed_tree, beta, B->tree.hashed_tree);
    }

    else if (tier == contiguous_tier)
    {
      A_matrix = spamm_chunk_get_matrix((*A)->tree.chunk);
      B_matrix = spamm_chunk_get_matrix(B->tree.chunk);

      for (i = 0; i < ipow(N_upper[0]-N_lower[0], 2); i++)
      {
        A_matrix[i] = beta*B_matrix[i];
        (*A)->norm = beta*B->norm;
        (*A)->norm2 = (*A)->norm*(*A)->norm;
      }
    }

    else
    {
      if ((*A)->tree.child == NULL)
      {
        (*A)->tree.child = calloc(ipow(2, number_dimensions), sizeof(struct spamm_recursive_node_t*));
      }

      spamm_recursive_copy(&(*A)->tree.child[0], beta, B->tree.child[0], number_dimensions, N_lower, N_upper, tier+1, contiguous_tier, use_linear_tree);
      spamm_recursive_copy(&(*A)->tree.child[1], beta, B->tree.child[1], number_dimensions, N_lower, N_upper, tier+1, contiguous_tier, use_linear_tree);
      spamm_recursive_copy(&(*A)->tree.child[2], beta, B->tree.child[2], number_dimensions, N_lower, N_upper, tier+1, contiguous_tier, use_linear_tree);
      spamm_recursive_copy(&(*A)->tree.child[3], beta, B->tree.child[3], number_dimensions, N_lower, N_upper, tier+1, contiguous_tier, use_linear_tree);
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

  int dim;

  unsigned int *N_lower;
  unsigned int *N_upper;

  assert(B != NULL);

  if (*A == NULL)
  {
    /* Create new matrix A. */
    *A = spamm_new(B->number_dimensions, B->N, B->contiguous_tier, B->N_block, B->use_linear_tree);
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

  if ((*A)->use_linear_tree != B->use_linear_tree)
  {
    SPAMM_FATAL("mismatch in use_linear_tree\n");
  }

  if ((*A)->contiguous_tier != B->contiguous_tier)
  {
    SPAMM_FATAL("mismatch in contiguous tier\n");
  }

  if (B->contiguous_tier == 0 && B->use_linear_tree)
  {
    spamm_hashed_copy(&(*A)->tree.hashed_tree, 1.0, B->tree.hashed_tree);
  }

  else
  {
    N_lower = calloc(B->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(B->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < B->number_dimensions; dim++)
    {
      N_upper[dim] = B->N_padded;
    }

    spamm_recursive_copy(&(*A)->tree.recursive_tree, 1.0, B->tree.recursive_tree,
        B->number_dimensions, N_lower, N_upper, 0, B->contiguous_tier, B->use_linear_tree);

    free(N_lower);
    free(N_upper);
  }
}
