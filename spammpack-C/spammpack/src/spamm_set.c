#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A)
{
  unsigned int tier, index, i_tier, j_tier, delta_index;
  unsigned int data_offset;
  unsigned int *index_key;
  GHashTable *node_hashtable;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    fprintf(stderr, "illegal index values for A_ij\n");
    exit(1);
  }

  /* Loop through tiers to construct the tree structure. */
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, tier));

    i_tier = i/delta_index;
    j_tier = j/delta_index;

    /* Construct linear index of the node on this tier. */
    index = spamm_index_2D(i_tier, j_tier);

    /* Get hash table at this tier. */
    node_hashtable = g_hash_table_lookup(A->tier_hashtable, &tier);

    if (tier < A->kernel_tier)
    {
      /* Check whether we already have a block at this tier. */
      if ((node = g_hash_table_lookup(node_hashtable, &index)) == NULL)
      {
        node = spamm_new_node(tier, index, spamm_index_3D_ik0(i_tier, j_tier), spamm_index_3D_0kj(i_tier, j_tier));

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
        data = spamm_new_block(tier, index, spamm_index_3D_ik0(i_tier, j_tier), spamm_index_3D_0kj(i_tier, j_tier));

        /* Allocate new tier index key. */
        index_key = (unsigned int*) malloc(sizeof(unsigned int));
        *index_key = index;

        /* Insert new data block into hashtable. */
        g_hash_table_insert(node_hashtable, index_key, data);
      }

      /* The data layout in this kernel tier dense matrix is broken into
       * blocks of 4x4 matrix blocks. The bocks are ordered in row-major
       * order, i.e.
       *
       * A_OFFSET_ij corresponds to the offset into the kernel tier matrix of
       * one basic 4x4 matrix.
       *
       * A_OFFSET_11  0*4*4
       * A_OFFSET_12  1*4*4
       * A_OFFSET_13  2*4*4
       * A_OFFSET_14  3*4*4
       * A_OFFSET_21  4*4*4
       * A_OFFSET_22  5*4*4
       * A_OFFSET_23  6*4*4
       * A_OFFSET_24  7*4*4
       * A_OFFSET_31  8*4*4
       * A_OFFSET_32  9*4*4
       * A_OFFSET_33 10*4*4
       * A_OFFSET_34 11*4*4
       * A_OFFSET_41 12*4*4
       * A_OFFSET_42 13*4*4
       * A_OFFSET_43 14*4*4
       * A_OFFSET_44 15*4*4
       */

      /* Calculate index on lowest tier, i.e. for basic 4x4 matrix blocks. */
      data_offset = spamm_index_kernel_block(i%delta_index, j%delta_index);

      data->block_dense[data_offset] = Aij;

      data->block_dense_dilated[4*data_offset+0] = Aij;
      data->block_dense_dilated[4*data_offset+1] = Aij;
      data->block_dense_dilated[4*data_offset+2] = Aij;
      data->block_dense_dilated[4*data_offset+3] = Aij;

      /* Update norms. */
    }
  }
}
