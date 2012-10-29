#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define NEW_NORM
#define SPAMM_SET_NO_ZERO

/** Recursively set a matrix element.
 *
 * @param number_dimensions The number of dimensions.
 * @param i An array of row/column indices.
 * @param Aij The value of the matrix element A(i,j).
 * @param N An array of matrix dimensions (unpadded).
 * @param N_lower An array of left-most column indices.
 * @param N_upper An array of right-most column indices.
 * @param tier The tier this node is on.
 * @param contiguous_tier The size of the contiguous submatrix block.
 * @param linear_tier The size of the submatrix that is stored in hashed format.
 * @param layout The layout of the matrix elements.
 * @param node The node.
 */
void
spamm_recursive_set (const unsigned int number_dimensions,
    const unsigned int *const i,
    const unsigned int *const N,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const unsigned int kernel_tier,
    const unsigned int N_block,
    const short use_linear_tree,
    const unsigned int depth,
    const float Aij,
    struct spamm_recursive_node_t **node)
{
  unsigned int *new_N_lower;
  unsigned int *new_N_upper;

  float *A;

  if (*node == NULL)
  {
    /* Allocate new node. */
    *node = spamm_recursive_new_node(tier, number_dimensions, contiguous_tier,
        N_block, use_linear_tree, N, N_lower, N_upper);
  }

  /* Update norm. */
  (*node)->norm2 += Aij*Aij;
  (*node)->norm   = sqrt((*node)->norm2);

  if (number_dimensions == 2 && tier == contiguous_tier && use_linear_tree)
  {
    if ((*node)->tree.hashed_tree == NULL)
    {
      (*node)->tree.hashed_tree = spamm_hashed_new(tier, kernel_tier, depth,
          N_lower[0], N_upper[0], N_lower[1], N_upper[1]);
    }
    spamm_hashed_set(i[0], i[1], Aij, (*node)->tree.hashed_tree);
  }

  else if (tier == contiguous_tier)
  {
    if ((*node)->tree.chunk == NULL)
    {
      (*node)->tree.chunk = spamm_new_chunk(number_dimensions, N_block, N,
          N_lower, N_upper);
    }

    /* Set matrix element. */
    A = spamm_chunk_get_matrix((*node)->tree.chunk);
    A[spamm_chunk_matrix_index(number_dimensions, N_block, N_lower, i)] = Aij;
  }

  else
  {
    /* Allocate children nodes. */
    if ((*node)->tree.child == NULL)
    {
      (*node)->tree.child = calloc(ipow(2, number_dimensions), sizeof(struct spamm_recursive_node_t*));
    }

    new_N_lower = calloc(number_dimensions, sizeof(unsigned int));
    new_N_upper = calloc(number_dimensions, sizeof(unsigned int));

    switch (number_dimensions)
    {
      case 1:
        if (i[0] < N_lower[0]+(N_upper[0]-N_lower[0])/2)
        {
          new_N_lower[0] = N_lower[0];
          new_N_upper[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2;

          spamm_recursive_set(number_dimensions, i, N, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, kernel_tier, N_block,
              use_linear_tree, depth, Aij, &((*node)->tree.child[0]));
        }

        else
        {
          new_N_lower[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2;
          new_N_upper[0] = N_upper[0];

          spamm_recursive_set(number_dimensions, i, N, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, kernel_tier, N_block,
              use_linear_tree, depth, Aij, &((*node)->tree.child[1]));
        }

        break;

      case 2:
        if (i[0] < N_lower[0]+(N_upper[0]-N_lower[0])/2 &&
            i[1] < N_lower[1]+(N_upper[1]-N_lower[1])/2)
        {
          new_N_lower[0] = N_lower[0];
          new_N_upper[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2;
          new_N_lower[1] = N_lower[1];
          new_N_upper[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2;

          spamm_recursive_set(number_dimensions, i, N, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, kernel_tier, N_block,
              use_linear_tree, depth, Aij, &((*node)->tree.child[0]));
        }

        else if (i[0] <  N_lower[0]+(N_upper[0]-N_lower[0])/2 &&
            i[1] >= N_lower[1]+(N_upper[1]-N_lower[1])/2)
        {
          new_N_lower[0] = N_lower[0];
          new_N_upper[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2;
          new_N_lower[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2;
          new_N_upper[1] = N_upper[1];

          spamm_recursive_set(number_dimensions, i, N, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, kernel_tier, N_block,
              use_linear_tree, depth, Aij, &((*node)->tree.child[1]));
        }

        else if (i[0] >= N_lower[0]+(N_upper[0]-N_lower[0])/2 &&
            i[1] <  N_lower[1]+(N_upper[1]-N_lower[1])/2)
        {
          new_N_lower[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2;
          new_N_upper[0] = N_upper[0];
          new_N_lower[1] = N_lower[1];
          new_N_upper[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2;

          spamm_recursive_set(number_dimensions, i, N, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, kernel_tier, N_block,
              use_linear_tree, depth, Aij, &((*node)->tree.child[2]));
        }

        else
        {
          new_N_lower[0] = N_lower[0]+(N_upper[0]-N_lower[0])/2;
          new_N_upper[0] = N_upper[0];
          new_N_lower[1] = N_lower[1]+(N_upper[1]-N_lower[1])/2;
          new_N_upper[1] = N_upper[1];

          spamm_recursive_set(number_dimensions, i, N, new_N_lower,
              new_N_upper, tier+1, contiguous_tier, kernel_tier, N_block,
              use_linear_tree, depth, Aij, &((*node)->tree.child[3]));
        }

        break;

      default:
        SPAMM_FATAL("not implemented\n");
    }

    free(new_N_lower);
    free(new_N_upper);
  }
}

/** Set an element in a matrix.
 *
 * @param i The row index.
 * @param j The column index.
 * @param Aij The value of the matrix element A(i,j).
 * @param A The matrix.
 */
void
spamm_hashed_set (const unsigned int i, const unsigned int j, const float Aij, struct spamm_hashed_t *A)
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
  struct spamm_hashed_node_t *node;
  struct spamm_hashed_data_t *data;

  float old_Aij = 0;

#ifdef NEW_NORM
  unsigned int next_tier;
  unsigned int child_index;
  struct spamm_hashed_node_t *child_node;
  struct spamm_hashed_data_t *child_data;
  struct spamm_hashtable_t *next_tier_hashtable;
#endif

  assert(A != NULL);

  if (i < A->M_lower || i >= A->M_upper)
  {
    SPAMM_FATAL("i (%u) out of bounding box [%u, %u)\n", i, A->M_lower, A->M_upper);
  }

  if (j < A->N_lower || j >= A->N_upper)
  {
    SPAMM_FATAL("j (%u) out of bounding box [%u, %u)\n", j, A->N_lower, A->N_upper);
  }

  /* In the trivial case, we simply return. */
#ifdef SPAMM_SET_NO_ZERO
  if (Aij == 0.0) { return; }
#endif

  /* Loop through tiers to construct the tree structure. */
  for (tier = A->tier; tier <= A->kernel_tier; tier++)
  {
    /* Calculate the size of the matrix block. */
    delta_index = (A->M_upper-A->M_lower)/(1 << (tier-A->tier));

    /* Calculate the matrix block indices. */
    i_tier = (i-A->M_lower)/delta_index;
    j_tier = (j-A->N_lower)/delta_index;

    /* Construct linear index of the node on this tier. */
    index = spamm_index_2D(i_tier, j_tier);

    /* Get hash table at this tier. */
    node_hashtable = A->tier_hashtable[tier-A->tier];

    if (tier < A->kernel_tier)
    {
      /* Check whether we already have a block at this tier. */
      if ((node = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        node = spamm_hashed_new_node(tier, index);
        spamm_hashtable_insert(node_hashtable, index, node);
      }
    }

    else
    {
      if ((data = spamm_hashtable_lookup(node_hashtable, index)) == NULL)
      {
        data = spamm_hashed_new_data(tier, index, A->layout);
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
      norm_offset = spamm_index_norm((i-A->M_lower)%delta_index/SPAMM_N_BLOCK, (j-A->N_lower)%delta_index/SPAMM_N_BLOCK);
      data_offset = spamm_index_kernel_block((i-A->M_lower)%delta_index, (j-A->N_lower)%delta_index, A->layout);
      data_offset_transpose = spamm_index_kernel_block_transpose((i-A->M_lower)%delta_index, (j-A->N_lower)%delta_index, A->layout);

#ifndef NEW_NORM
      /* For norm calculations, get original value of Aij. */
      old_Aij = data->block_dense[data_offset];
#endif

      /* Set new value. */
      data->block_dense[data_offset] = Aij;
      data->block_dense_store[data_offset] = Aij;
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
          old_Aij = data->block_dense[spamm_index_kernel_block_hierarchical(((i-A->M_lower)%delta_index)/SPAMM_N_BLOCK,
              ((j-A->N_lower)%delta_index)/SPAMM_N_BLOCK, i_basic, j_basic, A->layout)];
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

    }
  }

  /* Update norms back up the tree. Watch out for loop comparisons, tier is
   * unsigned and we better start looping top down. */
  for (tier = A->tier+1; tier <= A->kernel_tier; tier++)
  {
    reverse_tier = A->kernel_tier-tier+A->tier;

    delta_index = (A->M_upper-A->M_lower)/(1 << (reverse_tier-A->tier));

    i_tier = (i-A->M_lower)/delta_index;
    j_tier = (j-A->N_lower)/delta_index;

    /* Construct linear index of the node on this tier. */
    index = spamm_index_2D(i_tier, j_tier);

    /* Get hash table at this tier. */
    node_hashtable = A->tier_hashtable[reverse_tier-A->tier];

    /* Get the node. */
    node = spamm_hashtable_lookup(node_hashtable, index);

#ifdef NEW_NORM
    /* Get the tier hashtable for the next tier. */
    next_tier = node->tier+1;
    next_tier_hashtable = A->tier_hashtable[next_tier-A->tier];

    node->norm2 = 0.0;

    if (next_tier == A->kernel_tier)
    {
      for (i_child = 0; i_child < 2; i_child++) {
        for (j_child = 0; j_child < 2; j_child++)
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
      for (i_child = 0; i_child < 2; i_child++) {
        for (j_child = 0; j_child < 2; j_child++)
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

/** Set an element in a matrix.
 *
 * @param i The row/column index.
 * @param Aij The value of the matrix element A(i,j).
 * @param A The matrix.
 */
void
spamm_set (const unsigned int *const i, const float Aij, struct spamm_matrix_t *A)
{
  unsigned int *N_lower;
  unsigned int *N_upper;
  int dim;

  assert(A != NULL);

  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    if (i[dim] >= A->N[dim])
    {
      SPAMM_FATAL("i[%u] out of bounds (i[%u] = %i and N[%u] = %i)\n", dim, dim, i[dim], dim, A->N[dim]);
    }
  }

  /* Don't store zero. */
  if (Aij == 0.0) { return; }

  /* Store matrix element. */
  if (A->number_dimensions == 2 && A->contiguous_tier == 0 && A->use_linear_tree)
  {
    /* In case we only have a linear tree. */
    if (A->tree.hashed_tree == NULL)
    {
      A->tree.hashed_tree = spamm_hashed_new(0, A->kernel_tier, A->depth, 0, A->N_padded, 0, A->N_padded);
    }
    spamm_hashed_set(i[0], i[1], Aij, A->tree.hashed_tree);
  }

  else if (A->contiguous_tier == 0)
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    A->tree.chunk = spamm_new_chunk(A->number_dimensions, A->N_block, A->N, N_lower, N_upper);

    free(N_lower);
    free(N_upper);
  }

  else
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    spamm_recursive_set(A->number_dimensions, i, A->N, N_lower, N_upper, 0,
        A->contiguous_tier, A->kernel_tier, A->N_block, A->use_linear_tree,
        A->depth, Aij, &(A->tree.recursive_tree));

    free(N_lower);
    free(N_upper);
  }
}
