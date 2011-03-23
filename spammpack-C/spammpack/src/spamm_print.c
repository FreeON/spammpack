#include "spamm.h"
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
    const enum spamm_dense_type_t type, const float *A)
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
  struct spamm_node_t *node = value;

  printf("(node) tier %u: index_2D = %u, norm = %1.2e\n",
      node->tier, node->index_2D, node->norm);
}

void
spamm_print_data (unsigned int key, void *value, void *user_data)
{
  unsigned int i, j;
  struct spamm_data_t *data = value;

  printf("(node) tier %u: index_2D = %u, ", data->tier, data->index_2D);
  printf("node norm = %1.2e, ", data->node_norm);
  printf("norm = { ");
  for (i = 0; i < SPAMM_N_KERNEL_BLOCK; i++)
  {
    printf("{");
    for (j = 0; j < SPAMM_N_KERNEL_BLOCK; j++)
    {
      printf(" %1.2e", data->norm[i*SPAMM_N_KERNEL_BLOCK+j]);
    }
    printf(" }");
    if (i < SPAMM_N_KERNEL_BLOCK-1) { printf(", "); }
  }
  printf(" }, block_dense = ");
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
#ifdef SPAMM_USE_TRANSPOSE
  printf(" }, block_dense_transpose = ");
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
#endif
}

/** Print a SpAMM matrix, printing the tree structure.
 *
 * @param A The matrix.
 */
void
spamm_print_tree (const struct spamm_t *A)
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
spamm_print (const struct spamm_t *A)
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
