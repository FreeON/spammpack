#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Prune matrix tree by removing zero blocks.
 *
 * @param A The matrix.
 */
void
spamm_prune (struct spamm_hashed_t *A)
{
  spamm_error_fatal(__FILE__, __LINE__, "FIXME\n");
}

/** Expand a matrix tree to a full tree by adding zero blocks.
 *
 * @param A The matrix.
 */
void
spamm_expand (struct spamm_hashed_t *A)
{
  unsigned int i, j;
  unsigned int i_tier;
  unsigned int j_tier;
  unsigned int index;
  struct spamm_hashed_data_t *data;
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
        data = spamm_hashed_new_data(A->kernel_tier, index, A->layout);
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
spamm_construct_tree (struct spamm_hashed_t *A)
{
  unsigned int i;
  unsigned int tier;
  unsigned int next_tier;
  unsigned int reverse_tier;
  unsigned int parent_index;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_hashtable_t *next_tier_hashtable;
  struct spamm_list_t *tier_indices;
  struct spamm_hashed_data_t *data;
  struct spamm_hashed_node_t *node;
  struct spamm_hashed_node_t *parent_node;

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
          parent_node = spamm_hashed_new_node(next_tier, parent_index);
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
          parent_node = spamm_hashed_new_node(next_tier, parent_index);
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

/** Update the matrix norms.
 *
 * @param A The matrix.
 */
void
spamm_hashed_norm_update (struct spamm_hashed_t *A)
{
  int tier;
  int i_blocked, j_blocked;
  int i_basic, j_basic;
  int i_child, j_child;

  unsigned int i;
  unsigned int child_index;

  struct spamm_hashed_data_t *data;
  struct spamm_hashed_node_t *node;
  struct spamm_hashed_data_t *child_data;
  struct spamm_hashed_node_t *child_node;
  struct spamm_list_t *tier_index;

  assert(A != NULL);

  for (tier = A->kernel_tier; tier >= 0; tier--)
  {
    tier_index = spamm_hashtable_keys(A->tier_hashtable[tier]);

    if (tier == A->kernel_tier)
    {
      for (i = 0; i < spamm_list_length(tier_index); i++)
      {
        data = spamm_hashtable_lookup(A->tier_hashtable[tier], spamm_list_get_index(tier_index, i));

        /* Calculate norms on kernel blocks. */
        for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
          for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++)
          {
            data->norm2[spamm_index_norm(i_blocked, j_blocked)] = 0.0;
            for (i_basic = 0; i_basic < SPAMM_N_BLOCK; i_basic++) {
              for (j_basic = 0; j_basic < SPAMM_N_BLOCK; j_basic++)
              {
                data->norm2[spamm_index_norm(i_blocked, j_blocked)] +=
                  data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, data->layout)]
                  * data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, data->layout)];
              }
            }
            data->norm[spamm_index_norm(i_blocked, j_blocked)] = sqrt(data->norm2[spamm_index_norm(i_blocked, j_blocked)]);
          }
        }

        /* Calculate norms on kernel tier block. */
        data->node_norm2 = 0.0;
        for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
          for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++)
          {
            data->node_norm2 += data->norm2[spamm_index_norm(i_blocked, j_blocked)];
          }
        }
        data->node_norm = sqrt(data->node_norm2);
      }
    }

    else
    {
      for (i = 0; i < spamm_list_length(tier_index); i++)
      {
        node = spamm_hashtable_lookup(A->tier_hashtable[tier], spamm_list_get_index(tier_index, i));

        node->norm2 = 0.0;
        for (i_child = 0; i_child < 2; i_child++) {
          for (j_child = 0; j_child < 2; j_child++)
          {
            /* Construct index of child block. */
            child_index = (spamm_list_get_index(tier_index, i) << 2) | (i_child << 1) | j_child;

            if (tier+1 == A->kernel_tier)
            {
              child_data = spamm_hashtable_lookup(A->tier_hashtable[tier+1], child_index);

              if (child_data != NULL)
              {
                node->norm2 += child_data->node_norm2;
              }
            }

            else
            {
              child_node = spamm_hashtable_lookup(A->tier_hashtable[tier+1], child_index);

              if (child_node != NULL)
              {
                node->norm2 += child_node->norm2;
              }
            }
          }
        }
        node->norm = sqrt(node->norm2);
      }
    }
    spamm_list_delete(&tier_index);
  }
}
