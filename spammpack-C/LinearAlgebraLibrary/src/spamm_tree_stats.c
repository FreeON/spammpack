#include "spamm.h"
#include <assert.h>
#include <stdio.h>

void
spamm_node_stats (struct spamm_tree_stats_t *stats, const struct spamm_node_t *node)
{
  int i, j;

  stats->memory_tree += sizeof(struct spamm_node_t);
  if (node->child != NULL)
  {
    stats->memory_tree += node->M_child*node->N_child*sizeof(struct spamm_node_t*);
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        stats->number_nodes++;
        spamm_node_stats(stats, node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
      }
    }
  }

  if (node->block_dense != NULL)
  {
    stats->number_dense_blocks++;
    stats->memory_dense_blocks += node->M_block*node->N_block*sizeof(float_t);
  }
}

void
spamm_tree_stats (struct spamm_tree_stats_t *stats, const struct spamm_t *A)
{
  assert(stats != NULL);
  assert(A != NULL);

  stats->number_nodes = 0;
  stats->number_dense_blocks = 0;
  stats->memory_tree = sizeof(struct spamm_t);
  stats->memory_dense_blocks = 0;

  /* Recurse. */
  if (A->root != NULL)
  {
    stats->number_nodes++;
    spamm_node_stats(stats, A->root);
  }
}
