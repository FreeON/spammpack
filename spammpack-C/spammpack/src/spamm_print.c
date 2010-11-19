#include "spamm.h"
#include <assert.h>
#include <stdio.h>

void
spamm_print_node (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int *index = key;
  struct spamm_node_t *node = value;

  printf("tier %u: index = %u\n", node->tier, node->index);
}

void
spamm_print_data (gpointer key, gpointer value, gpointer user_data)
{
  unsigned int i, j;
  unsigned int *index = key;
  struct spamm_data_t *data = value;

  printf("tier %u: index = %u, block_dense = ", data->tier, data->index);
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
spamm_print (const struct spamm_t *A)
{
  unsigned int tier;
  GHashTable *node_hashtable;
  struct spamm_node_t *node;

  assert(A != NULL);

  printf("root node: M = %u, N = %u, N_padded = %u, depth = %u, kernel_tier = %u\n",
      A->M, A->N, A->N_padded, A->depth, A->kernel_tier);

  for (tier = 0; tier < A->kernel_tier; tier++)
  {
    if ((node_hashtable = g_hash_table_lookup(A->tier, &tier)) != NULL)
    {
      g_hash_table_foreach(node_hashtable, spamm_print_node, NULL);
    }
  }

  if ((node_hashtable = g_hash_table_lookup(A->tier, &A->kernel_tier)) != NULL)
  {
    g_hash_table_foreach(node_hashtable, spamm_print_data, NULL);
  }
}
