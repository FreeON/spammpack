#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Print a dense matrix.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param type The storage type of the dense matrix.
 * @param A The dense matrix.
 */
void
spamm_print_dense (const unsigned int M, const unsigned int N,
    const enum spamm_layout_t type, const float *A)
{
  unsigned int i, j;

  switch (type)
  {
    case column_major:
      for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++)
        {
          printf(" % 1.2e", A[spamm_index_column_major(i, j, M, N)]);
        }
        printf("\n");
      }
      break;

    case row_major:
      for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++)
        {
          printf(" % 1.2e", A[spamm_index_row_major(i, j, M, N)]);
        }
        printf("\n");
      }
      break;

    default:
      printf("error\n");
      exit(1);
      break;
  }
}

/** Print out a spamm_hashed_node_t object.
 *
 * @param key The key for the hashtable.
 * @param value The spamm_hashed_node_t object.
 * @param user_data Some additional data that can be passed into this function.
 */
void
spamm_print_node (unsigned int key, void *value, void *user_data)
{
  struct spamm_hashed_node_t *node = value;

  printf("(node) tier %u: index_2D = %u, norm = %1.2e\n",
      node->tier, node->index_2D, node->norm);
}

/** Print out a spamm_hashed_data_t object.
 *
 * @param key The key for the hashtable.
 * @param value The spamm_hashed_data_t object.
 * @param user_data Some additional data that can be passed into this function.
 */
void
spamm_print_data (unsigned int key, void *value, void *user_data)
{
  unsigned int i, j;
  struct spamm_hashed_data_t *data = value;

  printf("(node) tier %u: index_2D = %u, ", data->tier, data->index_2D);
  printf("node norm = %1.2e\n", data->node_norm);
  printf("norm = { ");
  for (i = 0; i < SPAMM_N_KERNEL_BLOCKED; i++)
  {
    printf("{");
    for (j = 0; j < SPAMM_N_KERNEL_BLOCKED; j++)
    {
      printf(" %1.2e", data->norm[i*SPAMM_N_KERNEL_BLOCKED+j]);
    }
    printf(" }");
    if (i < SPAMM_N_KERNEL_BLOCKED-1) { printf(", "); }
  }
  printf(" }, ");
  printf("block_dense = ");
  for (i = 0; i < SPAMM_N_KERNEL; i++)
  {
    printf("{");
    for (j = 0; j < SPAMM_N_KERNEL; j++)
    {
      printf(" %1.2e", data->block_dense[i*SPAMM_N_KERNEL+j]);
    }
    printf(" }");
    if (i < SPAMM_N_KERNEL-1) { printf(", "); }
  }
  printf(" }\n");
  printf("block_dense =\n");
  spamm_print_dense(SPAMM_N_KERNEL, SPAMM_N_KERNEL, row_major, data->block_dense);
  printf("block_dense_transpose = ");
  for (i = 0; i < SPAMM_N_KERNEL; i++)
  {
    printf("{");
    for (j = 0; j < SPAMM_N_KERNEL; j++)
    {
      printf(" %1.2e", data->block_dense_transpose[i*SPAMM_N_KERNEL+j]);
    }
    printf(" }");
    if (i < SPAMM_N_KERNEL-1) { printf(", "); }
  }
  printf(" }\n");
  printf("block_dense_transpose =\n");
  spamm_print_dense(SPAMM_N_KERNEL, SPAMM_N_KERNEL, row_major, data->block_dense_transpose);
}

/** Print a recursive tree node and all nodes underneath.
 *
 * @param node The tree node.
 */
void
spamm_recursive_print_node (const struct spamm_recursive_node_t *const node)
{
  unsigned int tier;
  unsigned int i;
  unsigned int dim;

  if (node == NULL) { return; }

  printf("node (tier = %u): N = {", node->tier);
  for (dim = 0; dim < node->number_dimensions; dim++)
  {
    printf(" [ %u --> %u ]", node->N_lower[dim], node->N_upper[dim]);
  }
  printf(" }, norm = %1.2e, ", node->norm);
  if (node->child == NULL)
  {
    printf("no children");
  }

  else
  {
    printf("{");
    for (i = 0; i < node->number_dimensions*node->number_dimensions; i++)
    {
      printf(" %p", node->child[i]);
    }
    printf(" }");
  }
  printf("\n");

  if (node->N_upper[0]-node->N_lower[0] == node->N_linear)
  {
    for (tier = node->tier; tier <= node->hashed_tree->kernel_tier; tier++)
    {
      if (tier < node->hashed_tree->kernel_tier)
      {
        spamm_hashtable_foreach(node->hashed_tree->tier_hashtable[tier-node->tier], spamm_print_node, NULL);
      }

      else
      {
        spamm_hashtable_foreach(node->hashed_tree->tier_hashtable[tier-node->tier], spamm_print_data, NULL);
      }
    }
  }

  else if (node->N_upper[0]-node->N_lower[0] == node->N_contiguous)
  {
    return;
  }

  else
  {
    switch (node->number_dimensions)
    {
      case 2:
        for (i = 0; i < 4; i++)
        {
          spamm_recursive_print_node(node->child[i]);
        }
        break;

      default:
        SPAMM_FATAL("not implemented\n");
        break;
    }
  }
}

/** Print a SpAMM matrix, printing the tree structure.
 *
 * @param A The matrix.
 */
void
spamm_print_tree (const struct spamm_matrix_t *const A)
{
  unsigned int dim;
  unsigned int tier;

  assert(A != NULL);

  printf("root node: ");
  printf("ndim = %u, ", A->number_dimensions);
  printf("N = {");
  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    printf(" %u", A->N[dim]);
  }
  printf(" }, ");
  printf("N_padded = %u, ", A->N_padded);
  printf("depth = %u, ", A->depth);
  printf("N_contiguous = %u, ", A->N_contiguous);
  printf("N_linear = %u, ", A->N_linear);
  printf("kernel_tier = %u, ", A->kernel_tier);
  printf("linear_tier = %u, ", A->linear_tier);
  printf("contiguous_tier = %u\n", A->contiguous_tier);

  if (A->recursive_tree != NULL)
  {
    spamm_recursive_print_node(A->recursive_tree);
  }

  else if (A->hashed_tree != NULL)
  {
    for (tier = 0; tier <= A->kernel_tier; tier++)
    {
      if (tier < A->kernel_tier)
      {
        spamm_hashtable_foreach(A->hashed_tree->tier_hashtable[tier], spamm_print_node, NULL);
      }

      else
      {
        spamm_hashtable_foreach(A->hashed_tree->tier_hashtable[tier], spamm_print_data, NULL);
      }
    }
  }
}

/** Print a SpAMM matrix in the same format a dense matrix is printed.
 *
 * @param A The matrix.
 */
void
spamm_print (const struct spamm_matrix_t *A)
{
  unsigned int *i;

  assert(A != NULL);

  i = calloc(A->number_dimensions, sizeof(unsigned int));
  switch(A->number_dimensions)
  {
    case 2:
      for (i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for (i[1] = 0; i[1] < A->N[1]; i[1]++)
        {
          printf(" % 1.2e", spamm_get(i, A));
        }
        printf("\n");
      }
      break;

    default:
      SPAMM_FATAL("not implemented\n");
  }
}
