#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

void
spamm_print_tree_node (const struct spamm_node_t *node)
{
  int i, j;

  spamm_print_node(node);
  if (node->child != NULL)
  {
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        spamm_print_tree_node(node->child[spamm_dense_index(i, j, node->N_child)]);
      }
    }
  }
}

void
spamm_print_tree (const struct spamm_t *A)
{
  printf("A: M = %i, N = %i, ", A->M, A->N);
  printf("M_padded = %i, N_padded = %i, ", A->M_padded, A->N_padded);
  printf("M_block = %i, N_block = %i, ", A->M_block, A->N_block);
  printf("M_child = %i, N_child = %i, ", A->M_child, A->N_child);
  printf("threshold = %7.1e, ", A->threshold);
  printf("root = %p\n", (void*) A->root);

  if (A->root != NULL)
  {
    spamm_print_tree_node(A->root);
  }
}
