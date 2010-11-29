#include "spamm.h"
#include <assert.h>
#include <stdio.h>

void
spamm_print_dense (const unsigned int M, const unsigned int N, const float *A)
{
  unsigned int i, j;

  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % 1.2f", A[spamm_index_row_major(i, j, M, N)]);
    }

    if (i < M-1)
    {
      printf("\n");
    }
  }
}

void
spamm_print_node (gpointer key, gpointer value, gpointer user_data)
{
  struct spamm_node_t *node = value;

  printf("(node) tier %u: index_2D = %u, index_3D_ik0 = %u, index_3D_0kj = %u\n",
      node->tier, node->index_2D, node->index_3D_ik0, node->index_3D_0kj);
}

void
spamm_print_data (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int i, j;
  struct spamm_data_t *data = value;

  printf("(node) tier %u: index_2D = %u, index_3D_ik0 = %u, index_3D_0kj = %u, ",
      data->tier, data->index_2D, data->index_3D_ik0, data->index_3D_0kj);
  printf("norm = { ");
  for (i = 0; i < SPAMM_N_KERNEL_BLOCK; i++)
  {
    printf("{");
    for (j = 0; j < SPAMM_N_KERNEL_BLOCK; j++)
    {
      printf(" %1.2f", data->norm[i*SPAMM_N_KERNEL_BLOCK+j]);
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
      printf(" %1.2f", data->block_dense[i*SPAMM_N_KERNEL+j]);
    }
    printf(" }");
    if (i < SPAMM_N_KERNEL-1) { printf(", "); }
  }
  printf(" }\n");
}

void
spamm_print_tier (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int *tier = key;
  unsigned int *kernel_tier = user_data;
  GHashTable *tier_hashtable = value;

  if (*tier < *kernel_tier)
  {
    g_hash_table_foreach(tier_hashtable, spamm_print_node, NULL);
  }

  else
  {
    g_hash_table_foreach(tier_hashtable, spamm_print_data, NULL);
  }
}

void
spamm_print (const struct spamm_t *A)
{
  assert(A != NULL);

  printf("root node: M = %u, N = %u, N_padded = %u, depth = %u, kernel_tier = %u\n",
      A->M, A->N, A->N_padded, A->depth, A->kernel_tier);

  g_hash_table_foreach(A->tier_hashtable, spamm_print_tier, (gpointer) &A->kernel_tier);
}
