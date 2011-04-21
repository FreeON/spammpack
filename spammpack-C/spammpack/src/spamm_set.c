#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NEW_NORM
#define SPAMM_SET_NO_ZERO

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
  unsigned int i_blocked;
  unsigned int j_blocked;
  unsigned int i_basic;
  unsigned int j_basic;
  unsigned int i_child;
  unsigned int j_child;
  unsigned int delta_index;
  unsigned int norm_offset;
  unsigned int data_offset;
  unsigned int data_offset_transpose;
  struct spamm_hashtable_t *node_hashtable;
  struct spamm_node_t *node;
  struct spamm_data_t *data;

  float old_Aij = 0;

  float norm_A11, norm_A12, norm_A21, norm_A22;

#ifdef NEW_NORM
  unsigned int next_tier;
  unsigned int child_index;
  struct spamm_node_t *child_node;
  struct spamm_data_t *child_data;
  struct spamm_hashtable_t *next_tier_hashtable;
#endif

  assert(A != NULL);
  assert(i < A->M);
  assert(j < A->N);

  if (i >= A->M || j >= A->N)
  {
    fprintf(stderr, "illegal index values for A_ij\n");
    exit(1);
  }

  /* In the trivial case, we simply return. */
#ifdef SPAMM_SET_NO_ZERO
  if (Aij == 0.0) { return; }
#endif

  /* Loop through tiers to construct the tree structure. */
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    /* Calculate the size of the matrix block. */
    delta_index = (unsigned int) floor(A->N_padded/pow(SPAMM_N_CHILD, tier));

    /* Calculate the matrix block indices. */
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
        node = spamm_new_node(tier, index);
        spamm_hashtable_insert(node_hashtable, index, node);
      }
    }

    else
    {
      if ((data = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        data = spamm_new_block(tier, index, A->layout);
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

      /* Calculate offsets into the norm and the matrix data. */
      norm_offset = spamm_index_norm(i%delta_index/SPAMM_N_BLOCK, j%delta_index/SPAMM_N_BLOCK);
      data_offset = spamm_index_kernel_block(i%delta_index, j%delta_index, A->layout);
      data_offset_transpose = spamm_index_kernel_block_transpose(i%delta_index, j%delta_index, A->layout);

#ifndef NEW_NORM
      /* For norm calculations, get original value of Aij. */
      old_Aij = data->block_dense[data_offset];
#endif

      /* Set new value. */
      data->block_dense[data_offset] = Aij;
      data->block_dense_transpose[data_offset_transpose] = Aij;

      data->block_dense_dilated[4*data_offset+0] = Aij;
      data->block_dense_dilated[4*data_offset+1] = Aij;
      data->block_dense_dilated[4*data_offset+2] = Aij;
      data->block_dense_dilated[4*data_offset+3] = Aij;

      /* Update norms. */
#ifdef NEW_NORM
      data->norm2[norm_offset] = 0.0;
      for (i_basic = 0; i_basic < SPAMM_N_BLOCK; i_basic++) {
        for (j_basic = 0; j_basic < SPAMM_N_BLOCK; j_basic++)
        {
          old_Aij = data->block_dense[spamm_index_kernel_block_hierarchical((i%delta_index)/SPAMM_N_BLOCK, (j%delta_index)/SPAMM_N_BLOCK, i_basic, j_basic, A->layout)];
          data->norm2[norm_offset] += old_Aij*old_Aij;
        }
      }
      data->norm[norm_offset] = sqrt(data->norm2[norm_offset]);

      data->node_norm2 = 0.0;
      for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
        for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++)
        {
          data->node_norm2 += data->norm2[spamm_index_norm(i_blocked, j_blocked)];
        }
      }
      data->node_norm = sqrt(data->node_norm2);
#else
      data->norm2[norm_offset] += Aij*Aij-old_Aij*old_Aij;
      data->norm[norm_offset] = sqrt(data->norm2[norm_offset]);

      data->node_norm2 += Aij*Aij-old_Aij*old_Aij;
      data->node_norm = sqrt(data->node_norm2);
#endif

      /* Update upper tier norms. */
      norm_A11 = sqrt(
          data->norm2[spamm_index_norm(0, 0)]+
          data->norm2[spamm_index_norm(0, 1)]+
          data->norm2[spamm_index_norm(1, 0)]+
          data->norm2[spamm_index_norm(1, 1)]);
      norm_A12 = sqrt(
          data->norm2[spamm_index_norm(0, 2)]+
          data->norm2[spamm_index_norm(0, 3)]+
          data->norm2[spamm_index_norm(1, 2)]+
          data->norm2[spamm_index_norm(1, 3)]);
      norm_A21 = sqrt(
          data->norm2[spamm_index_norm(2, 0)]+
          data->norm2[spamm_index_norm(2, 1)]+
          data->norm2[spamm_index_norm(3, 0)]+
          data->norm2[spamm_index_norm(3, 1)]);
      norm_A22 = sqrt(
          data->norm2[spamm_index_norm(2, 2)]+
          data->norm2[spamm_index_norm(2, 3)]+
          data->norm2[spamm_index_norm(3, 2)]+
          data->norm2[spamm_index_norm(3, 3)]);

      data->norm_upper[0] = norm_A11;
      data->norm_upper[1] = norm_A12;
      data->norm_upper[2] = norm_A11;
      data->norm_upper[3] = norm_A12;
      data->norm_upper[4] = norm_A21;
      data->norm_upper[5] = norm_A22;
      data->norm_upper[6] = norm_A21;
      data->norm_upper[7] = norm_A22;

      data->norm_upper_transpose[0] = norm_A11;
      data->norm_upper_transpose[1] = norm_A21;
      data->norm_upper_transpose[2] = norm_A12;
      data->norm_upper_transpose[3] = norm_A22;
      data->norm_upper_transpose[4] = norm_A11;
      data->norm_upper_transpose[5] = norm_A21;
      data->norm_upper_transpose[6] = norm_A12;
      data->norm_upper_transpose[7] = norm_A22;
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

#ifdef NEW_NORM
    /* Get the tier hashtable for the next tier. */
    next_tier = node->tier+1;
    next_tier_hashtable = A->tier_hashtable[next_tier];

    node->norm2 = 0.0;

    if (next_tier == A->kernel_tier)
    {
      for (i_child = 0; i_child < SPAMM_N_CHILD; i_child++) {
        for (j_child = 0; j_child < SPAMM_N_CHILD; j_child++)
        {
          /* Construct index of child block. */
          child_index = (index << 2) | (i_child << 1) | j_child;

          /* Get child node. */
          child_data = spamm_hashtable_lookup(next_tier_hashtable, child_index);

          if (child_data != NULL)
          {
            node->norm2 += child_data->node_norm2;
          }
        }
      }
    }

    else
    {
      for (i_child = 0; i_child < SPAMM_N_CHILD; i_child++) {
        for (j_child = 0; j_child < SPAMM_N_CHILD; j_child++)
        {
          /* Construct index of child block. */
          child_index = (index << 2) | (i_child << 1) | j_child;

          /* Get child node. */
          child_node = spamm_hashtable_lookup(next_tier_hashtable, child_index);

          if (child_node != NULL)
          {
            node->norm2 += child_node->norm2;
          }
        }
      }
    }
    node->norm = sqrt(node->norm2);
#else
    /* Update norms. */
    node->norm2 += Aij*Aij-old_Aij*old_Aij;
    node->norm = sqrt(node->norm2);
#endif
  }
}
