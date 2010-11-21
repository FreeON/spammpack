#include "spamm.h"
#include <assert.h>
#include <stdio.h>

void
spamm_print_node (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int *index = key;
  struct spamm_node_t *node = value;

  printf("(node) tier %u: index = %u\n", node->tier, node->index);
}

void
spamm_print_data (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int i, j;
  unsigned int *index = key;
  struct spamm_data_t *data = value;

  printf("(data) tier %u: index = %u, block_dense = ", data->tier, data->index);
  for (i = 0; i < SPAMM_N_KERNEL; i++) {
    printf("{");
    for (j = 0; j < SPAMM_N_KERNEL; j++)
    {
      printf(" %1.2f", data->block_dense[i][j]);
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
