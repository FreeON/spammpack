#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

/** Print detailed information of a tree node.
 *
 * This function prints detailed information of a tree node. If prints out the
 * values of all fields in struct spamm_node_t.
 */
void
spamm_print_node (const struct spamm_node_t *node)
{
  int i, j;
  int nonzero;

  if (node == NULL)
  {
    printf("node %p\n", (void*) node);
  }

  else
  {
    printf("node %p: ", (void*) node);
    printf("M_lower = %i, M_upper = %i, ", node->M_lower, node->M_upper);
    printf("N_lower = %i, N_upper = %i, ", node->N_lower, node->N_upper);
    printf("M_block = %i, N_block = %i, ", node->M_block, node->N_block);
    printf("M_child = %i, N_child = %i, ", node->M_child, node->N_child);
    printf("threshold = %7.1e, ", node->threshold);
    printf("index = %u, ", node->index);
    printf("child = %p, ", (void*) node->child);
    if (node->child != NULL)
    {
      printf("{ ");
      for (i = 0; i < node->M_child; ++i) {
        for (j = 0; j < node->N_child; ++j)
        {
          printf("%p", (void*) node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
          if (j < node->N_child-1) { printf(", "); }
        }
        if (i < node->M_child-1) { printf(", "); }
      }
      printf(" }, ");
    }
    printf("block_dense = %p", (void*) node->block_dense);
    if (node->block_dense != NULL)
    {
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
      printf(", sparsity = %1.1f%%", (float_t) nonzero / (float_t) (node->M_block*node->N_block));
      printf(", { ");
      for (i = 0; i < node->M_block; ++i) {
        for (j = 0; j < node->N_block; ++j)
        {
          printf("%f", node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)]);
          if (j < node->N_block-1) { printf(", "); }
        }
        if (i < node->M_block-1) { printf(", "); }
      }
      printf(" }\n");
    }

    else { printf("\n"); }
  }
}
