/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 */
void
spamm_hashed_add (const float alpha,
    struct spamm_hashed_t **A,
    const float beta,
    struct spamm_hashed_t **B)
{
  struct spamm_list_t *A_keys;
  struct spamm_list_t *B_keys;

  struct spamm_hashed_node_t *node;
  struct spamm_hashed_data_t *A_data;
  struct spamm_hashed_data_t *B_data;

  struct spamm_hashed_node_t *child_node;
  struct spamm_hashed_data_t *child_data;

  unsigned int A_index;
  unsigned int B_index;
  unsigned int child_index;

  int i, j, k;
  int i_blocked, j_blocked;
  int i_basic, j_basic;
  int i_child, j_child;

  int tier;
  int next_tier;

  float Aij;

  if (*A == NULL && *B == NULL) { return; }

  if (*A != NULL && *B != NULL)
  {
    if ((*A)->layout != (*B)->layout)
    {
      SPAMM_FATAL("inconsistent layout in matrices\n");
    }

    if ((*A)->M_upper-(*A)->M_lower != (*B)->M_upper-(*B)->M_lower)
    {
      SPAMM_FATAL("mismatch of number of rows\n");
    }

    if ((*A)->N_upper-(*A)->N_lower != (*B)->N_upper-(*B)->N_lower)
    {
      SPAMM_FATAL("mismatch of number of columns\n");
    }

    /* Go to lowest tier and start adding submatrix blocks. */
    A_keys = spamm_hashtable_keys((*A)->tier_hashtable[(*A)->kernel_tier-(*A)->tier]);
    B_keys = spamm_hashtable_keys((*B)->tier_hashtable[(*A)->kernel_tier-(*A)->tier]);

    /* Sort the keys. */
    spamm_list_sort_index(A_keys, spamm_list_compare_int);
    spamm_list_sort_index(B_keys, spamm_list_compare_int);

    for (i = 0, j = 0; ; )
    {
      if (i < spamm_list_length(A_keys))
      {
        A_index = spamm_list_get_index(A_keys, i);
      }

      else
      {
        A_index = UINT_MAX;
      }

      if (j < spamm_list_length(B_keys))
      {
        B_index = spamm_list_get_index(B_keys, j);
      }

      else
      {
        B_index = UINT_MAX;
      }

      if (i >= spamm_list_length(A_keys) && j >= spamm_list_length(B_keys))
      {
        /* Done. */
        break;
      }

      if (A_index < B_index)
      {
        A_data = spamm_hashtable_lookup((*A)->tier_hashtable[(*A)->kernel_tier-(*A)->tier], A_index);
        for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
          for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++) {
            for (i_basic = 0; i_basic < SPAMM_N_BLOCK; i_basic++) {
              for (j_basic = 0; j_basic < SPAMM_N_BLOCK; j_basic++)
              {
                Aij = alpha*A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)];

                A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_transpose[spamm_index_kernel_block_transpose_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_store[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+0] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+1] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+2] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+3] = Aij;
              }
            }
          }
        }
        A_data->node_norm2 *= alpha*alpha;
        A_data->node_norm *= alpha;

        i++;
      }

      else if (B_index < A_index)
      {
        B_data = spamm_hashtable_lookup((*B)->tier_hashtable[(*B)->kernel_tier-(*B)->tier], spamm_list_get_index(B_keys, j));

        /* Create new node. */
        A_data = spamm_hashed_new_data((*A)->kernel_tier-(*A)->tier, B_index, (*A)->layout);
        A_data->node_norm2 = 0;
        for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
          for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++) {
            for (i_basic = 0; i_basic < SPAMM_N_BLOCK; i_basic++) {
              for (j_basic = 0; j_basic < SPAMM_N_BLOCK; j_basic++)
              {
                Aij = beta*B_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)];

                A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_transpose[spamm_index_kernel_block_transpose_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_store[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+0] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+1] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+2] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+3] = Aij;

                A_data->node_norm2 += Aij*Aij;
              }
            }
          }
        }
        A_data->node_norm = sqrt(A_data->node_norm2);

        j++;

        /* Create new block in A and store it. */
        spamm_hashtable_insert((*A)->tier_hashtable[(*A)->kernel_tier-(*A)->tier], B_index, A_data);
      }

      else if (A_index == B_index)
      {
        A_data = spamm_hashtable_lookup((*A)->tier_hashtable[(*A)->kernel_tier-(*A)->tier], spamm_list_get_index(A_keys, i));
        B_data = spamm_hashtable_lookup((*B)->tier_hashtable[(*B)->kernel_tier-(*B)->tier], spamm_list_get_index(B_keys, j));

        /* Add block data. */
        A_data->node_norm2 = 0;
        for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
          for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++) {
            for (i_basic = 0; i_basic < SPAMM_N_BLOCK; i_basic++) {
              for (j_basic = 0; j_basic < SPAMM_N_BLOCK; j_basic++)
              {
                Aij = alpha*A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)]
                  +beta*B_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)];

                A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_transpose[spamm_index_kernel_block_transpose_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_store[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+0] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+1] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+2] = Aij;
                A_data->block_dense_dilated[4*spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, i_basic, j_basic, (*A)->layout)+3] = Aij;

                A_data->node_norm2 += Aij*Aij;
              }
            }
          }
        }
        A_data->node_norm = sqrt(A_data->node_norm2);

        i++;
        j++;
      }

      else
      {
        SPAMM_FATAL("I should not be here\n");
      }

      /* Update block norms. */
      for (i_blocked = 0; i_blocked < SPAMM_N_KERNEL_BLOCKED; i_blocked++) {
        for (j_blocked = 0; j_blocked < SPAMM_N_KERNEL_BLOCKED; j_blocked++)
        {
          A_data->norm2[spamm_index_norm(i_blocked, j_blocked)] = 0;
          for (k = 0; k < SPAMM_N_BLOCK*SPAMM_N_BLOCK; k++)
          {
            A_data->norm2[spamm_index_norm(i_blocked, j_blocked)] +=
              A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, 0, 0, (*A)->layout)+k]
              *A_data->block_dense[spamm_index_kernel_block_hierarchical(i_blocked, j_blocked, 0, 0, (*A)->layout)+k];
          }
          A_data->norm[spamm_index_norm(i_blocked, j_blocked)] = sqrt(A_data->norm2[spamm_index_norm(i_blocked, j_blocked)]);
        }
      }
    }
    spamm_list_delete(&A_keys);
    spamm_list_delete(&B_keys);

    /* Fix norms of upper tiers. */
    for (tier = (*A)->kernel_tier-(*A)->tier-1; tier >= 0; tier--)
    {
      next_tier = tier+1;

      /* Construct missing nodes one tier up. */
      A_keys = spamm_hashtable_keys((*A)->tier_hashtable[next_tier]);
      for (i = 0; i < spamm_list_length(A_keys); i++)
      {
        A_index = spamm_list_get_index(A_keys, i) >> 2;
        if (spamm_hashtable_lookup((*A)->tier_hashtable[tier], A_index) == NULL)
        {
          spamm_hashtable_insert((*A)->tier_hashtable[tier], A_index, spamm_hashed_new_node(tier, A_index));
        }
      }
      spamm_list_delete(&A_keys);

      /* Update norms. */
      A_keys = spamm_hashtable_keys((*A)->tier_hashtable[tier]);
      for (i = 0; i < spamm_list_length(A_keys); i++)
      {
        A_index = spamm_list_get_index(A_keys, i);
        node = spamm_hashtable_lookup((*A)->tier_hashtable[tier], A_index);
        node->norm2 = 0;

        if (next_tier == (*A)->kernel_tier-(*A)->tier)
        {
          for (i_child = 0; i_child < 2; i_child++) {
            for (j_child = 0; j_child < 2; j_child++)
            {
              /* Construct index of child block. */
              child_index = (A_index << 2) | (i_child << 1) | j_child;

              /* Get child node. */
              child_data = spamm_hashtable_lookup((*A)->tier_hashtable[next_tier], child_index);

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
              child_index = (A_index << 2) | (i_child << 1) | j_child;

              /* Get child node. */
              child_node = spamm_hashtable_lookup((*A)->tier_hashtable[next_tier], child_index);

              if (child_node != NULL)
              {
                node->norm2 += child_node->norm2;
              }
            }
          }
        }
        node->norm = sqrt(node->norm2);
      }
      spamm_list_delete(&A_keys);
    }
  }

  else if (*A == NULL && *B != NULL)
  {
    SPAMM_FATAL("[FIXME]\n");
  }

  else if (*A != NULL && *B == NULL)
  {
    SPAMM_FATAL("[FIXME]\n");
  }
}

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 * @param number_dimensions The number of dimensions.
 */
void
spamm_recursive_add (const float alpha,
    struct spamm_recursive_node_t **A,
    const float beta,
    struct spamm_recursive_node_t **B,
    const unsigned int number_dimensions,
    const unsigned int *const N_lower,
    const unsigned int *const N_upper,
    const unsigned int tier,
    const unsigned int contiguous_tier,
    const short use_linear_tree)
{
  unsigned int i;

  float *A_matrix;
  float *B_matrix;

  /* There is nothing to do here. */
  if ((*A) == NULL && (*B) == NULL) { return; }

  if ((*A) != NULL && (*B) != NULL)
  {
    if (tier == contiguous_tier && use_linear_tree)
    {
      spamm_hashed_add(alpha, &(*A)->tree.hashed_tree, beta, &(*B)->tree.hashed_tree);
    }

    else if (tier == contiguous_tier)
    {
      if ((*A)->tree.chunk == NULL)
      {
        SPAMM_FATAL("??\n");
      }

      if ((*B)->tree.chunk == NULL)
      {
        SPAMM_FATAL("??\n");
      }

      /* Braindead add. */
      A_matrix = spamm_chunk_get_matrix((*A)->tree.chunk);
      B_matrix = spamm_chunk_get_matrix((*B)->tree.chunk);

      switch (number_dimensions)
      {
        case 2:
          for (i = 0; i < ipow(N_upper[0]-N_lower[0], 2); i++)
          {
            A_matrix[i] = alpha*A_matrix[i]+beta*B_matrix[i];
          }
          break;

        default:
          SPAMM_FATAL("not implemented\n");
      }
    }

    else
    {
      if ((*A)->tree.child == NULL)
      {
        SPAMM_FATAL("??\n");
      }

      if ((*B)->tree.child == NULL)
      {
        SPAMM_FATAL("??\n");
      }

      /* Recurse. */
      for (i = 0; i < ipow(2, number_dimensions); i++)
      {
        spamm_recursive_add(alpha, &(*A)->tree.child[i], beta,
            &(*B)->tree.child[i], number_dimensions, N_lower, N_upper,
            tier+1, contiguous_tier, use_linear_tree);
      }
    }
  }

  else if ((*A) == NULL && (*B) != NULL)
  {
    /* Copy B node to A. */
    spamm_recursive_copy(&(*A), beta, (*B), number_dimensions, N_lower, N_upper, tier, contiguous_tier, use_linear_tree);
  }

  else if ((*A) != NULL && (*B) == NULL)
  {
    /* Multiply A by alpha. */
    spamm_recursive_multiply_scalar(alpha, *A, number_dimensions, tier, contiguous_tier, use_linear_tree);
  }
}

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 */
void
spamm_add (const float alpha,
    struct spamm_matrix_t *A,
    const float beta,
    struct spamm_matrix_t *B)
{
  int dim;

  unsigned int *N_lower;
  unsigned int *N_upper;

  assert(A != NULL);
  assert(B != NULL);

  if (A->number_dimensions != B->number_dimensions)
  {
    SPAMM_FATAL("mismatch of number of dimensions\n");
  }

  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    if (A->N[dim] != B->N[dim])
    {
      SPAMM_FATAL("mismatch of number of rows/columns\n");
    }
  }

  if (A->contiguous_tier == 0 && A->use_linear_tree)
  {
    if (A->tree.hashed_tree != NULL || B->tree.hashed_tree != NULL)
    {
      spamm_hashed_add(alpha, &A->tree.hashed_tree, beta, &B->tree.hashed_tree);
    }
  }

  else if (A->tree.recursive_tree != NULL || B->tree.recursive_tree != NULL)
  {
    N_lower = calloc(A->number_dimensions, sizeof(unsigned int));
    N_upper = calloc(A->number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < A->number_dimensions; dim++)
    {
      N_upper[dim] = A->N_padded;
    }

    spamm_recursive_add(alpha, &A->tree.recursive_tree, beta,
        &B->tree.recursive_tree, A->number_dimensions, N_lower, N_upper, 0,
        A->contiguous_tier, A->use_linear_tree);

    free(N_lower);
    free(N_upper);
  }

  else
  {
    /* Matrices are empty. */
  }
}
