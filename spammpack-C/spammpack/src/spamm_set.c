#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A)
{
  unsigned int tier, i_tier, j_tier, index;
  GHashTable *node_hashtable;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    printf("illegal index values for A_ij\n");
    exit(1);
  }

  /* Loop through tiers to construct the tree structure. */
  for (tier = 0; tier < A->kernel_tier; tier++)
  {
    /* Construct linear index of the node on this tier. */
    i_tier = 0;
    j_tier = 0;
    index = 0;

    /* Check whether we already have a hash table at this tier. */
    if ((node_hashtable = g_hash_table_lookup(A->tier, &tier)) == NULL)
    {
      /* Create new hash table for this tier. */
      node_hashtable = g_hash_table_new(g_int_hash, g_int_equal);
      g_hash_table_insert(A->tier, &tier, node_hashtable);
    }

    /* Check whether we already have a block at this tier. */
    if ((node = g_hash_table_lookup(node_hashtable, &index)) == NULL)
    {
      node = spamm_new_node();
      g_hash_table_insert(node_hashtable, &index, node);
    }
  }

  /* Put the matrix element into the right place. */
  if ((node_hashtable = g_hash_table_lookup(A->tier, &A->kernel_tier)) == NULL)
  {
    node_hashtable = g_hash_table_new(g_int_hash, g_int_equal);
    g_hash_table_insert(A->tier, &A->kernel_tier, node_hashtable);
  }

  if ((data = g_hash_table_lookup(node_hashtable, &index)) == NULL)
  {
    data = spamm_new_block(A->kernel_tier, index);
    g_hash_table_insert(node_hashtable, &index, data);
  }
}
