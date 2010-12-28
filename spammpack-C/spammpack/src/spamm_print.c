#include "spamm.h"
#include <assert.h>
#include <stdio.h>

/** Print a dense matrix.
 *
 * @param M The number of rows.
 * @param N The number of columns.
 * @param A The dense matrix.
 */
void
spamm_print_dense (const unsigned int M, const unsigned int N, const float *A)
{
  unsigned int i, j;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % 1.2e", A[spamm_index_row_major(i, j, M, N)]);
    }
    printf("\n");
  }
}

void
spamm_print_node (unsigned int key, void *value, void *user_data)
{
  struct spamm_node_t *node = value;

  printf("(node) tier %u: index_2D = %u, index_3D_ik0 = %u, index_3D_0kj = %u, norm = %1.2e\n",
      node->tier, node->index_2D, node->index_3D_ik0, node->index_3D_0kj, node->norm);
}

void
spamm_print_data (unsigned int key, void *value, void *user_data)
{
  unsigned int i, j;
  struct spamm_data_t *data = value;

  printf("(node) tier %u: index_2D = %u, index_3D_ik0 = %u, index_3D_0kj = %u, ",
      data->tier, data->index_2D, data->index_3D_ik0, data->index_3D_0kj);
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
}

/** Print a SpAMM matrix.
 *
 * @param A The matrix.
 */
void
spamm_print (const struct spamm_t *A)
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
