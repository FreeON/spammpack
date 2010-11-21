#include "spamm.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A)
{
  unsigned int tier, index, i_tier, j_tier, delta_index;
  unsigned int *tier_key;
  unsigned int *index_key;
  GHashTable *node_hashtable;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    printf("illegal index values for A_ij\n");
    exit(1);
  }

  //printf("setting A(%u,%u)\n", i, j);

  /* Loop through tiers to construct the tree structure. */
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, tier));

    i_tier = (unsigned int) floor(i/delta_index);
    j_tier = (unsigned int) floor(j/delta_index);

    //printf("tier %u: delta = %u, (%u,%u) --> (%u,%u)\n", tier, delta_index, i, j, i_tier, j_tier);

    /* Construct linear index of the node on this tier. */
    index = spamm_index_2D(i_tier, j_tier);

    /* Get hash table at this tier. */
    node_hashtable = g_hash_table_lookup(A->tier_hashtable, &tier);

    if (tier < A->kernel_tier)
    {
      /* Check whether we already have a block at this tier. */
      if ((node = g_hash_table_lookup(node_hashtable, &index)) == NULL)
      {
        //printf("creating new node for tier %u with index %u\n", tier, index);
        node = spamm_new_node(tier, index);

        /* Allocate new linear index key. */
        index_key = (unsigned int*) malloc(sizeof(unsigned int));
        *index_key = index;

        /* Insert new linear index key into hashtable. */
        g_hash_table_insert(node_hashtable, index_key, node);
      }
    }

    else
    {
      if ((data = g_hash_table_lookup(node_hashtable, &index)) == NULL)
      {
        //printf("creating new data for tier %u with index %u\n", tier, index);
        data = spamm_new_block(tier, index);

        /* Allocate new tier index key. */
        index_key = (unsigned int*) malloc(sizeof(unsigned int));
        *index_key = index;

        /* Insert new data block into hashtable. */
        g_hash_table_insert(node_hashtable, index_key, data);
      }

      //printf("setting data at (%u,%u) and index %u\n", i-i_tier*delta_index, j-j_tier*delta_index, index);
      data->block_dense[i-i_tier*delta_index][j-j_tier*delta_index] = Aij;
    }

    //printf("tier %u: size(node_hashtable) = %u\n", tier, g_hash_table_size(node_hashtable));
  }
}
