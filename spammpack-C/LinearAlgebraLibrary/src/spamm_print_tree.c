#include "spamm.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** \private Print out a tree node and its children nodes.
 *
 * This is the recursive part of spamm_print_tree().
 *
 * @param node The node to print.
 * @param width The width of the binary representation of the linear index.
 */
void
spamm_print_tree_node (const struct spamm_node_t *node, const int width)
{
  int i, j;

  assert(node != NULL);

  spamm_print_node(node);
  for (i = 0; i < SPAMM_N_CHILD; ++i) {
    for (j = 0; j < SPAMM_N_CHILD; ++j)
    {
      if (node->child[i][j] != NULL)
      {
        spamm_print_tree_node(node->child[i][j], width);
      }
    }
  }
}

/** Print out the tree structure of a SpAMM matrix.
 *
 * @param A The matrix.
 */
void
spamm_print_tree (const struct spamm_t *A)
{
  assert(A != NULL);

  printf("A: M = %u, N = %u, ", A->M, A->N);
  printf("M_padded = %u, N_padded = %u, ", A->M_padded, A->N_padded);
  printf("M_block = %u, N_block = %u, ", SPAMM_N_BLOCK, SPAMM_N_BLOCK);
  printf("M_child = %u, N_child = %u, ", SPAMM_N_CHILD, SPAMM_N_CHILD);
  printf("depth = %u, ", A->tree_depth);
  printf("root = %p\n", (void*) A->root);

  if (A->root != NULL)
  {
    spamm_print_tree_node(A->root, A->tree_depth*2);
  }
}
