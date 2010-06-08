#include "spamm.h"
#include <assert.h>
#include <stdio.h>

/** \private Statistics of a node.
 *
 * This is the recursive part, use spamm_tree_stats() instead.
 *
 * @param stats The spamm_tree_stats_t structure with the result.
 * @param node The node to check.
 */
void
spamm_node_stats (struct spamm_tree_stats_t *stats, const struct spamm_node_t *node)
{
  int i, j;
  int nonzero;

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

    /* Calculate sparsity of dense block. */
    nonzero = 0;
    for (i = 0; i < node->M_block; ++i) {
      for (j = 0; j < node->N_block; ++j)
      {
        if (node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)] != 0.0)
        {
          nonzero++;
        }
      }
    }
    stats->average_sparsity += 1 - (float_t) nonzero / (float_t) (node->M_block*node->N_block);
  }
}

/** Get statistics about the tree.
 *
 * This function returns a struct spamm_tree_stats_t with information about
 * the tree.
 *
 * @param stats The spamm_tree_stats_t structure with the result.
 * @param A The matrix to check.
 */
void
spamm_tree_stats (struct spamm_tree_stats_t *stats, const struct spamm_t *A)
{
  assert(stats != NULL);
  assert(A != NULL);

  stats->number_nodes = 0;
  stats->number_dense_blocks = 0;
  stats->memory_tree = sizeof(struct spamm_t);
  stats->memory_dense_blocks = 0;
  stats->average_sparsity = 0;

  /* Recurse. */
  if (A->root != NULL)
  {
    stats->number_nodes++;
    spamm_node_stats(stats, A->root);
    stats->average_sparsity /= (float_t) stats->number_dense_blocks;
  }
}
