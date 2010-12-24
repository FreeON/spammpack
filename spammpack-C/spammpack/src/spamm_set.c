#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NEW_NORM

/** Set an element in a matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param Aij The value of the matrix element A(i,j).
 * @param A The matrix.
 */
void
spamm_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_t *A)
{
  unsigned int tier;
  unsigned int reverse_tier;
  unsigned int index;
  unsigned int i_tier;
  unsigned int j_tier;
  unsigned int i_block;
  unsigned int j_block;
  unsigned int delta_index;
  unsigned int norm_offset;
  unsigned int data_offset;
  unsigned int *index_key;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  float old_Aij = 0;

  assert(A != NULL);

  if (i >= A->M || j >= A->N)
  {
    fprintf(stderr, "illegal index values for A_ij\n");
    exit(1);
  }

  /* In the trivial case, we simply return. */
  if (Aij == 0.0) { return; }

  /* Loop through tiers to construct the tree structure. */
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, tier));

    i_tier = i/delta_index;
    j_tier = j/delta_index;

    /* Construct linear index of the node on this tier. */
    index = spamm_index_2D(i_tier, j_tier);

    /* Get hash table at this tier. */
    node_hashtable = A->tier_hashtable[tier];

    if (tier < A->kernel_tier)
    {
      /* Check whether we already have a block at this tier. */
      if ((node = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        node = spamm_new_node(tier, index, spamm_index_3D_ik0(i_tier, j_tier), spamm_index_3D_0kj(i_tier, j_tier));
        spamm_hashtable_insert(node_hashtable, index, node);
      }
    }

    else
    {
      if ((data = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        data = spamm_new_block(tier, index, spamm_index_3D_ik0(i_tier, j_tier), spamm_index_3D_0kj(i_tier, j_tier));
        spamm_hashtable_insert(node_hashtable, index, data);
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
      norm_offset = spamm_index_row_major((i%delta_index)/SPAMM_N_KERNEL_BLOCK, (j%delta_index)/SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK);
      data_offset = spamm_index_kernel_block(i%delta_index, j%delta_index);

      /* For norm calculations, get original value of Aij. */
      old_Aij = data->block_dense[data_offset];

      /* Set new value. */
      data->block_dense[data_offset] = Aij;

      data->block_dense_dilated[4*data_offset+0] = Aij;
      data->block_dense_dilated[4*data_offset+1] = Aij;
      data->block_dense_dilated[4*data_offset+2] = Aij;
      data->block_dense_dilated[4*data_offset+3] = Aij;

      /* Update norms. */
#ifdef NEW_NORM
      data->norm2[norm_offset] = 0.0;
      for (i_block = 0; i_block < SPAMM_N_BLOCK; i_block++) {
        for (j_block = 0; j_block < SPAMM_N_BLOCK; j_block++)
        {
          old_Aij = data->block_dense[SPAMM_N_BLOCK*SPAMM_N_BLOCK
            *spamm_index_row_major(i%delta_index/SPAMM_N_KERNEL_BLOCK,
                j%delta_index/SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)
            +spamm_index_row_major(i_block, j_block, SPAMM_N_BLOCK, SPAMM_N_BLOCK)];
          data->norm2[norm_offset] += old_Aij*old_Aij;
        }
      }
      data->norm[norm_offset] = sqrt(data->norm2[norm_offset]);

      data->node_norm2 = 0.0;
      for (i_block = 0; i_block < SPAMM_N_KERNEL_BLOCK; i_block++) {
        for (j_block = 0; j_block < SPAMM_N_KERNEL_BLOCK; j_block++)
        {
          data->node_norm2 += data->norm2[spamm_index_row_major(i_block, j_block, SPAMM_N_KERNEL_BLOCK, SPAMM_N_KERNEL_BLOCK)];
        }
      }
      data->node_norm = sqrt(data->node_norm2);
#else
      data->norm2[norm_offset] += Aij*Aij-old_Aij*old_Aij;
      data->norm[norm_offset] = sqrt(data->norm2[norm_offset]);

      data->node_norm2 += Aij*Aij-old_Aij*old_Aij;
      data->node_norm = sqrt(data->node_norm2);
#endif
    }
  }

  /* Update norms back up the tree. Watch out for loop comparisons, tier is
   * unsigned and we better start looping top down. */
  for (tier = 1; tier <= A->kernel_tier; tier++)
  {
    reverse_tier = A->kernel_tier-tier;

    delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, reverse_tier));

    i_tier = i/delta_index;
    j_tier = j/delta_index;

    /* Construct linear index of the node on this tier. */
    index = spamm_index_2D(i_tier, j_tier);

    /* Get hash table at this tier. */
    node_hashtable = A->tier_hashtable[reverse_tier];

    /* Get the node. */
    node = spamm_hashtable_lookup(node_hashtable, index);

    /* Update norms. */
    node->norm2 += Aij*Aij-old_Aij*old_Aij;
    node->norm = sqrt(node->norm2);
  }
}
