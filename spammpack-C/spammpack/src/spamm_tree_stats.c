#include "spamm.h"
#include <assert.h>
#include <math.h>
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
  unsigned int kernel_block_N;

  assert(node != NULL);

  stats->memory_tree += sizeof(struct spamm_node_t);
  stats->memory_tree += SPAMM_N_CHILD*SPAMM_N_CHILD*sizeof(struct spamm_node_t*);
  for (i = 0; i < SPAMM_N_CHILD; ++i) {
    for (j = 0; j < SPAMM_N_CHILD; ++j)
    {
      if (node->child[i][j] != NULL)
      {
        stats->number_nodes++;
        spamm_node_stats(stats, node->child[i][j]);
      }
    }
  }

  if (node->tier == node->kernel_tier)
  {
    kernel_block_N = pow(SPAMM_N_CHILD, node->tree_depth-node->kernel_tier)*SPAMM_N_BLOCK;

    stats->number_dense_blocks++;
    stats->memory_dense_blocks += kernel_block_N*kernel_block_N*sizeof(floating_point_t);

    /* Calculate sparsity of dense block. */
    nonzero = 0;
    for (i = 0; i < kernel_block_N; ++i) {
      for (j = 0; j < kernel_block_N; ++j)
      {
        if (node->block_dense[spamm_dense_index(i, j, kernel_block_N, kernel_block_N)] != 0.0)
        {
          nonzero++;
        }
      }
    }
    stats->average_sparsity += 1 - (floating_point_t) nonzero / (floating_point_t) (kernel_block_N*kernel_block_N);
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
    stats->average_sparsity /= (floating_point_t) stats->number_dense_blocks;
  }

  stats->memory_total = stats->memory_tree+stats->memory_dense_blocks;
}
