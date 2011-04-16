#include "spamm.h"
#include <math.h>

/** Prune matrix tree by removing zero blocks.
 *
 * @param A The matrix.
 */
void
spamm_prune (struct spamm_t *A)
{
}

/** Expand a matrix tree to a full tree by adding zero blocks.
 *
 * @param A The matrix.
 */
void
spamm_expand (struct spamm_t *A)
{
  unsigned int i, j;
  unsigned int i_tier;
  unsigned int j_tier;
  unsigned int index;
  struct spamm_data_t *data;
  struct spamm_hashtable_t *node_hashtable;

  /* Get the kernel tier hashtable. */
  node_hashtable = A->tier_hashtable[A->kernel_tier];

  /* Create new nodes at the kernel tier. */
  for (i = 0; i < A->M; i += SPAMM_N_KERNEL) {
    for (j = 0; j < A->N; j += SPAMM_N_KERNEL)
    {
      /* Calculate the matrix block indices. */
      i_tier = i/SPAMM_N_KERNEL;
      j_tier = j/SPAMM_N_KERNEL;

      /* Construct linear index of the node on this tier. */
      index = spamm_index_2D(i_tier, j_tier);

      /* Create new block in case one doesn't exist already. */
      if ((data = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        data = spamm_new_block(A->kernel_tier, index, A->layout);
        spamm_hashtable_insert(node_hashtable, index, data);
      }
    }
  }

  /* Construct the rest of the tree. */
  spamm_construct_tree(A);
}

/** Construct the full tree from blocks at the kernel tier level.
 *
 * @param A The matrix.
 */
void
spamm_construct_tree (struct spamm_t *A)
{
  unsigned int i;
  unsigned int tier;
  unsigned int next_tier;
  unsigned int reverse_tier;
  unsigned int parent_index;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_hashtable_t *next_tier_hashtable;
  struct spamm_list_t *tier_indices;
  struct spamm_data_t *data;
  struct spamm_node_t *node;
  struct spamm_node_t *parent_node;

  if (A->kernel_tier == 0) { return; }

  /* Create tree above kernel tier. */
  for (tier = 0; tier <= A->kernel_tier-1; tier++)
  {
    reverse_tier = A->kernel_tier-tier;
    node_hashtable = A->tier_hashtable[reverse_tier];
    next_tier = reverse_tier-1;
    next_tier_hashtable = A->tier_hashtable[next_tier];

    tier_indices = spamm_hashtable_keys(node_hashtable);

    if (reverse_tier == A->kernel_tier)
    {
      for (i = 0; i < spamm_list_length(tier_indices); i++)
      {
        /* Get block. */
        data = spamm_hashtable_lookup(node_hashtable, spamm_list_get_index(tier_indices, i));

        /* Construct index of parent node. */
        parent_index = data->index_2D >> 2;

        /* Get parent node. */
        parent_node = spamm_hashtable_lookup(next_tier_hashtable, parent_index);

        if (parent_node == NULL)
        {
          parent_node = spamm_new_node(next_tier, parent_index);
          spamm_hashtable_insert(next_tier_hashtable, parent_index, parent_node);
        }

        parent_node->norm2 += data->node_norm2;
      }
    }

    else
    {
      for (i = 0; i < spamm_list_length(tier_indices); i++)
      {
        /* Get block. */
        node = spamm_hashtable_lookup(node_hashtable, spamm_list_get_index(tier_indices, i));

        /* Construct index of parent node. */
        parent_index = node->index_2D >> 2;

        /* Get parent node. */
        parent_node = spamm_hashtable_lookup(next_tier_hashtable, parent_index);

        if (parent_node == NULL)
        {
          parent_node = spamm_new_node(next_tier, parent_index);
          spamm_hashtable_insert(next_tier_hashtable, parent_index, parent_node);
        }

        parent_node->norm2 += node->norm2;
      }
    }
    spamm_list_delete(&tier_indices);

    /* Update norms. */
    tier_indices = spamm_hashtable_keys(next_tier_hashtable);
    for (i = 0; i < spamm_list_length(tier_indices); i++)
    {
      /* Get block. */
      node = spamm_hashtable_lookup(next_tier_hashtable, spamm_list_get_index(tier_indices, i));
      node->norm = sqrt(node->norm2);
    }
    spamm_list_delete(&tier_indices);
  }
}
