#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

#define RELATIVE_TOLERANCE 1.0e-7

struct spamm_check_user_data_t
{
  unsigned int tier;
  const struct spamm_t *A;
};

void
spamm_check_verify_norm (gpointer key, gpointer value, gpointer user_data)
{
  short i, j;
  short i_block, j_block;
  float norm2 = 0.0;

  float Aij;

  unsigned int child_index;
  struct spamm_node_t *child_node = NULL;
  struct spamm_data_t *child_data = NULL;

  unsigned int *index = key;
  struct spamm_node_t *node = NULL;
  struct spamm_data_t *data = NULL;
  struct spamm_check_user_data_t *user = user_data;

  unsigned int next_tier;
  GHashTable *next_tier_hashtable = NULL;

  /* Load correct value. */
  if (user->tier == user->A->kernel_tier)
  {
    data = value;

    /* Check norms on kernel blocks. */
    for (i = 0; i < SPAMM_N_KERNEL_BLOCK; i++) {
      for (j = 0; j < SPAMM_N_KERNEL_BLOCK; j++)
      {
        norm2 = 0.0;
        for (i_block = 0; i_block < SPAMM_N_BLOCK; i_block++) {
          for (j_block = 0; j_block < SPAMM_N_BLOCK; j_block++)
          {
            Aij = data->block_dense[SPAMM_N_BLOCK*SPAMM_N_BLOCK*spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)
              +spamm_index_row_major(i_block, j_block, SPAMM_N_BLOCK, SPAMM_N_BLOCK)];
            norm2 += Aij*Aij;
          }
        }

        if (fabs(norm2-data->norm2[spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)]) > RELATIVE_TOLERANCE*norm2)
        {
          printf("tier %u, index %u, block (%u,%u): incorrect norm value, found %e, should be %e, |diff| = %e, fixing...\n",
              data->tier, data->index_2D, i, j,
              data->norm[spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)],
              sqrt(norm2),
              fabs(data->norm[spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)]-sqrt(norm2)));

          data->norm2[spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)] = norm2;
          data->norm[spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)] = sqrt(norm2);
        }
      }
    }

    /* Check norms on kernel tier block. */
    norm2 = 0.0;
    for (i = 0; i < SPAMM_N_KERNEL; i++) {
      for (j = 0; j < SPAMM_N_KERNEL; j++)
      {
        Aij = data->block_dense[spamm_index_row_major(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)];
        norm2 += Aij*Aij;
      }
    }

    if (fabs(norm2-data->node_norm2) > RELATIVE_TOLERANCE*norm2)
    {
      printf("tier %u, index %u: incorrect node norm value, found %e, should be %e, |diff| = %e, fixing...\n",
          data->tier, data->index_2D, data->node_norm, sqrt(norm2), fabs(data->node_norm-sqrt(norm2)));

      data->node_norm2 = norm2;
      data->node_norm = sqrt(norm2);
    }
  }

  else
  {
    node = value;

    /* Get the tier hashtable for the next tier. */
    next_tier = user->tier+1;
    next_tier_hashtable = g_hash_table_lookup(user->A->tier_hashtable, &next_tier);

    if (next_tier == user->A->kernel_tier)
    {
      norm2 = 0.0;
      for (i = 0; i < SPAMM_N_CHILD; i++) {
        for (j = 0; j < SPAMM_N_CHILD; j++)
        {
          /* Construct index of child block. */
          child_index = ((*index) << 2) | (i << 1) | j;

          child_data = g_hash_table_lookup(next_tier_hashtable, &child_index);

          norm2 += child_data->node_norm2;
        }
      }
    }

    else
    {
      norm2 = 0.0;
      for (i = 0; i < SPAMM_N_CHILD; i++) {
        for (j = 0; j < SPAMM_N_CHILD; j++)
        {
          /* Construct index of child block. */
          child_index = ((*index) << 2) | (i << 1) | j;

          child_node = g_hash_table_lookup(next_tier_hashtable, &child_index);

          norm2 += child_node->norm2;
        }
      }
    }

    if (fabs(norm2-node->norm2) > RELATIVE_TOLERANCE*norm2)
    {
      printf("tier %u, index %u: incorrect norm value, found %e, should be %e, |diff| = %e, fixing...\n",
          node->tier, node->index_2D, node->norm, sqrt(norm2), fabs(node->norm-sqrt(norm2)));

      node->norm2 = norm2;
      node->norm = sqrt(norm2);
    }
  }
}

/** Check the internal consistency of a matrix.
 *
 * @param A The matrix to check
 *
 * @return The following error codes are returned:
 *   - SPAMM_OK - The matrix is consistent.
 *   - SPAMM_ERROR - Something is not consistent.
 */
int
spamm_check (const struct spamm_t *A)
{
  int result = SPAMM_OK;
  unsigned int depth;
  unsigned int N_padded;
  unsigned int tier;
  unsigned int reverse_tier;
  float x_M, x_N, x;
  struct spamm_check_user_data_t user_data;
  GHashTable *hashtable;

  assert(A != NULL);

  /* Calculate the padding and depth of matrix based on values stored in M and
   * N.
   */
  x_M = (log(A->M) > log(SPAMM_N_BLOCK) ? log(A->M) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);
  x_N = (log(A->N) > log(SPAMM_N_BLOCK) ? log(A->N) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (depth >= 1 && ((int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, depth-1)) >= A->M && (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, depth-1)) >= A->N))
  {
    depth--;
  }

  if (A->depth != depth)
  {
    printf("depth incorrect, should be %u, but is %u\n", depth, A->depth);
    return SPAMM_ERROR;
  }

  N_padded = (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, depth));

  if (A->N_padded != N_padded)
  {
    printf("padded matrix dimensions incorrect, should be %u, but is %u\n", N_padded, A->N_padded);
    return SPAMM_ERROR;
  }

  if (A->kernel_tier != depth-SPAMM_KERNEL_DEPTH)
  {
    printf("kernel tier incorrect, should be %u, but is %u\n", depth-SPAMM_KERNEL_DEPTH, A->kernel_tier);
    return SPAMM_ERROR;
  }

  /* Check whether there are tier hashtables for every tier. */
  if (A->tier_hashtable == NULL)
  {
    printf("no tier hashtable\n");
    return SPAMM_ERROR;
  }

  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    if ((hashtable = g_hash_table_lookup(A->tier_hashtable, &tier)) == NULL)
    {
      printf("missing tier hashtable for tier %u\n", tier);
      return SPAMM_ERROR;
    }
  }

  /* Check each node. */
  for (tier = A->kernel_tier+1; tier >= 1; tier--)
  {
    /* Get tier hashtable. */
    reverse_tier = tier-1;
    hashtable = g_hash_table_lookup(A->tier_hashtable, &reverse_tier);

    /* Verify consistency of 2D and 3D linear indices. */

    /* Verify norms. */
    user_data.tier = reverse_tier;
    user_data.A = A;
    g_hash_table_foreach(hashtable, spamm_check_verify_norm, &user_data);
  }

  return result;
}
