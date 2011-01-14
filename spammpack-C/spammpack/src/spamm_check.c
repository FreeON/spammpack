#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define RELATIVE_TOLERANCE 1.0e-7

/** @private The user data object that is passed to to the hash table
 * iterator.
 */
struct spamm_check_user_data_t
{
  /** The result. */
  int result;

  /** The tier this hash table is on. */
  unsigned int tier;

  /** The spamm_t object the hash table was taken from. */
  const struct spamm_t *A;
};

/** Verify tree structure.
 *
 * @param index The linear index of this node.
 * @param value The matrix tree node. This can have either type spamm_node_t
 * or spamm_data_t depending on the tier.
 * @param user_data The user data of type spamm_check_user_data_t.
 */
void
spamm_check_tree_structure (unsigned int index, void *value, void *user_data)
{
}

/** Verify linear indices.
 *
 * @param index The linear index of this node.
 * @param value The matrix tree node. This can have either type spamm_node_t
 * or spamm_data_t depending on the tier.
 * @param user_data The user data of type spamm_check_user_data_t.
 */
void
spamm_check_linear_index (unsigned int index, void *value, void *user_data)
{
}

/** Verify the norm of a node.
 *
 * This function is called in the hash table iterator for each tier of the
 * matrix tree.
 *
 * @param index The linear index of this node.
 * @param value The matrix tree node. This can have either type spamm_node_t
 * or spamm_data_t depending on the tier.
 * @param user_data The user data of type spamm_check_user_data_t.
 */
void
spamm_check_norm (unsigned int index, void *value, void *user_data)
{
  short i, j;
  short i_block, j_block;
  float norm2 = 0.0;

  float Aij;

  unsigned int child_index;
  struct spamm_node_t *child_node = NULL;
  struct spamm_data_t *child_data = NULL;

  struct spamm_node_t *node = NULL;
  struct spamm_data_t *data = NULL;
  struct spamm_check_user_data_t *user = user_data;

  unsigned int next_tier;
  struct spamm_hashtable_t *next_tier_hashtable = NULL;

  unsigned int norm_offset;

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
            Aij = data->block_dense[SPAMM_N_BLOCK*SPAMM_N_BLOCK
              *spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)
              +spamm_index_row_major(i_block, j_block, SPAMM_N_BLOCK, SPAMM_N_BLOCK)];
            norm2 += Aij*Aij;
          }
        }

        norm_offset = spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK);
        if (fabs(norm2-data->norm2[norm_offset]) > RELATIVE_TOLERANCE*norm2 ||
            fabs(sqrt(norm2)-data->norm[norm_offset]) > RELATIVE_TOLERANCE*sqrt(norm2))
        {
          printf("tier %u, index %u, block (%u,%u): incorrect block norm value, found %e, should be %e, |diff| = %e, rel. diff = %e, fixing...\n",
              data->tier, data->index_2D, i, j,
              data->norm[norm_offset],
              sqrt(norm2),
              fabs(data->norm[norm_offset]-sqrt(norm2)),
              (norm2 != 0.0 ? fabs(data->norm[norm_offset]-sqrt(norm2))/sqrt(norm2) : 0));

          data->norm2[norm_offset] = norm2;
          data->norm[norm_offset] = sqrt(norm2);
          user->result = SPAMM_ERROR;
        }
      }
    }

    /* Check norms on kernel tier block. */
    norm2 = 0.0;
    for (i = 0; i < SPAMM_N_KERNEL_BLOCK; i++) {
      for (j = 0; j < SPAMM_N_KERNEL_BLOCK; j++)
      {
        norm2 += data->norm2[spamm_index_row_major(i, j, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)];
      }
    }

    if (fabs(norm2-data->node_norm2) > RELATIVE_TOLERANCE*norm2 ||
        fabs(sqrt(norm2)-data->node_norm) > RELATIVE_TOLERANCE*sqrt(norm2))
    {
      printf("tier %u, index %u: incorrect node norm value, found %e = sqrt(%e), should be %e, |diff| = %e, rel. diff = %e, fixing...\n",
          data->tier, data->index_2D,
          data->node_norm,
          data->node_norm2,
          sqrt(norm2),
          fabs(data->node_norm-sqrt(norm2)),
          (norm2 != 0.0 ? fabs(data->node_norm-sqrt(norm2))/sqrt(norm2) : 0));

      data->node_norm2 = norm2;
      data->node_norm = sqrt(norm2);
      user->result = SPAMM_ERROR;
    }
  }

  else
  {
    node = value;

    /* Get the tier hashtable for the next tier. */
    next_tier = user->tier+1;
    next_tier_hashtable = user->A->tier_hashtable[next_tier];

    norm2 = 0.0;

    if (next_tier == user->A->kernel_tier)
    {
      for (i = 0; i < SPAMM_N_CHILD; i++) {
        for (j = 0; j < SPAMM_N_CHILD; j++)
        {
          /* Construct index of child block. */
          child_index = (index << 2) | (i << 1) | j;

          /* Get child node. */
          child_data = spamm_hashtable_lookup(next_tier_hashtable, child_index);

          if (child_data != NULL)
          {
            norm2 += child_data->node_norm2;
          }
        }
      }
    }

    else
    {
      for (i = 0; i < SPAMM_N_CHILD; i++) {
        for (j = 0; j < SPAMM_N_CHILD; j++)
        {
          /* Construct index of child block. */
          child_index = (index << 2) | (i << 1) | j;

          /* Get child node. */
          child_node = spamm_hashtable_lookup(next_tier_hashtable, child_index);

          if (child_node != NULL)
          {
            norm2 += child_node->norm2;
          }
        }
      }
    }

    if (fabs(norm2-node->norm2) > RELATIVE_TOLERANCE*norm2 ||
        fabs(sqrt(norm2)-node->norm) > RELATIVE_TOLERANCE*sqrt(norm2))
    {
      printf("tier %u, index %u: incorrect norm value, found %e = sqrt(%e), should be %e, |diff| = %e, rel. diff = %e, fixing...\n",
          node->tier, node->index_2D,
          node->norm,
          node->norm2,
          sqrt(norm2),
          fabs(node->norm-sqrt(norm2)),
          (norm2 != 0.0 ? fabs(node->norm-sqrt(norm2))/sqrt(norm2) : 0));

      node->norm2 = norm2;
      node->norm = sqrt(norm2);
      user->result = SPAMM_ERROR;
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
  unsigned int depth;
  unsigned int N_padded;
  unsigned int tier;
  unsigned int reverse_tier;
  float x_M, x_N, x;
  struct spamm_check_user_data_t user_data;
  struct spamm_hashtable_t *hashtable;

  assert(A != NULL);

  user_data.result = SPAMM_OK;

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
    if ((hashtable = A->tier_hashtable[tier]) == NULL)
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
    hashtable = A->tier_hashtable[reverse_tier];

    user_data.tier = reverse_tier;
    user_data.A = A;

    /* Verify tree structure. */
    spamm_hashtable_foreach(hashtable, spamm_check_tree_structure, &user_data);

    /* Verify consistency of linear indices. */
    spamm_hashtable_foreach(hashtable, spamm_check_linear_index, &user_data);

    /* Verify norms. */
    spamm_hashtable_foreach(hashtable, spamm_check_norm, &user_data);
  }

  return user_data.result;
}
