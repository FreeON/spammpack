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

void
spamm_print_node (unsigned int key, void *value, void *user_data)
{
  struct spamm_hashed_node_t *node = value;

  printf("(node) tier %u: index_2D = %u, norm = %1.2e\n",
      node->tier, node->index_2D, node->norm);
}

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
  printf("norm_upper = { ");
  for (i = 0; i < 8; i++)
  {
    printf(" %1.2e", data->norm_upper[i]);
  }
  printf(" }, ");
  printf("norm_upper_transpose = { ");
  for (i = 0; i < 8; i++)
  {
    printf(" %1.2e", data->norm_upper_transpose[i]);
  }
  printf(" }\n");
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

/** Print a SpAMM matrix, printing the tree structure.
 *
 * @param A The matrix.
 */
void
spamm_print_tree (const struct spamm_hashed_t *A)
{
  unsigned int tier;

  assert(A != NULL);

  printf("root node: M = %u, N = %u, N_padded = %u, depth = %u, kernel_tier = %u\n",
      A->M, A->N, A->N_padded, A->depth, A->kernel_tier);

  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    if (tier < A->kernel_tier)
    {
      spamm_hashtable_foreach(A->tier_hashtable[tier], spamm_print_node, NULL);
    }

    else
    {
      spamm_hashtable_foreach(A->tier_hashtable[tier], spamm_print_data, NULL);
    }
  }
}

/** Print a SpAMM matrix in the same format a dense matrix is printed.
 *
 * @param A The matrix.
 */
void
spamm_print (const struct spamm_hashed_t *A)
{
  unsigned int i;
  unsigned int j;

  assert(A != NULL);

  for (i = 0; i < A->M; i++) {
    for (j = 0; j < A->N; j++)
    {
      printf(" % 1.2e", spamm_get(i, j, A));
    }
    printf("\n");
  }
}
