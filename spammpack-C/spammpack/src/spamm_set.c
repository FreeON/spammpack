#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A)
{
  unsigned int tier, i_tier, j_tier, index, delta_index;
  GHashTable *node_hashtable;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    printf("illegal index values for A_ij\n");
    exit(1);
  }

  printf("setting A(%u,%u)\n", i, j);

  /* Loop through tiers to construct the tree structure. */
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    /* Construct linear index of the node on this tier. */
    delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, tier));

    i_tier = (unsigned int) floor(i/delta_index);
    j_tier = (unsigned int) floor(j/delta_index);

    printf("tier %u: delta = %u, (%u,%u) --> (%u,%u)\n", tier, delta_index, i, j, i_tier, j_tier);

    index = spamm_index_2D(i_tier, j_tier);

    /* Check whether we already have a hash table at this tier. */
    node_hashtable = g_hash_table_lookup(A->tier, &tier);
    if (node_hashtable == NULL)
    {
      /* Create new hash table for this tier. */
      printf("creating new hashtable for tier %u\n", tier);
      node_hashtable = g_hash_table_new(g_int_hash, spamm_hash_uint_equal);
      g_hash_table_insert(A->tier, &tier, node_hashtable);
    }

    if (tier < A->kernel_tier)
    {
      /* Check whether we already have a block at this tier. */
      node = g_hash_table_lookup(node_hashtable, &index);
      if (node == NULL)
      {
        printf("creating new node for tier %u with index %u\n", tier, index);
        node = spamm_new_node(tier, index);
        g_hash_table_insert(node_hashtable, &index, node);
      }
    }

    else
    {
      data = g_hash_table_lookup(node_hashtable, &index);
      if (data == NULL)
      {
        printf("creating new data for tier %u with index %u\n", tier, index);
        data = spamm_new_block(tier, index);
        g_hash_table_insert(node_hashtable, &index, data);
      }

      printf("setting data at index %u\n", index);
    }

    printf("tier %u: size(node_hashtable) = %u\n", tier, g_hash_table_size(node_hashtable));
  }
}
