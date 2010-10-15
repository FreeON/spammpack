#include "spamm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** Print detailed information of a tree node.
 *
 * This function prints detailed information of a tree node. If prints out the
 * values of all fields in struct spamm_node_t.
 *
 * @param node The matrix node.
 */
void
spamm_print_node (const struct spamm_node_t *node)
{
  int i, j;
  int nonzero;
  unsigned int kernel_block_M, kernel_block_N;
  char *binary_string;
  char header[1000];
  char empty_header[1000];

  if (node == NULL)
  {
    printf("node %p\n", (void*) node);
  }

  else
  {
    sprintf(header, "node %p: ", (void*) node);
    for (i = 0; i < strlen(header); ++i)
    {
      empty_header[i] = ' ';
    }
    empty_header[strlen(header)] = '\0';
    printf("%s", header);
    printf("tier = %i, ", node->tier);
    printf("depth = %i, ", node->tree_depth);
    printf("linear = %i, ", node->linear_tier);
    printf("A(%u:%u,%u:%u), ", node->M_lower+1, node->M_upper, node->N_lower+1, node->N_upper);
    printf("M_lower = %i, M_upper = %i, ", node->M_lower, node->M_upper);
    printf("N_lower = %i, N_upper = %i, ", node->N_lower, node->N_upper);
    printf("\n");
    printf("%s", empty_header);
    binary_string = (char*) malloc(sizeof(char)*(node->tier*2+2));
    spamm_int_to_binary(node->index, node->tier*2, binary_string);
    printf("index = 0b%s (%u), ", binary_string, node->index);
    free(binary_string);
    printf("child = %p, ", (void*) node->child);
    if (node->child != NULL)
    {
      printf("(1:%u,1:%u) { ", SPAMM_M_CHILD, SPAMM_N_CHILD);
      for (i = 0; i < SPAMM_M_CHILD; ++i) {
        printf("(%u,1:%u) { ", i+1, SPAMM_N_CHILD);
        for (j = 0; j < SPAMM_N_CHILD; ++j)
        {
          printf("%p", (void*) node->child[i][j]);
          if (j < SPAMM_N_CHILD-1) { printf(", "); }
        }
        printf(" }");
        if (i < SPAMM_M_CHILD-1) { printf(", "); }
      }
      printf(" }\n");
      printf("%s", empty_header);
    }
    printf("linear = %p", (void*) node->linear_quadtree);
    if (node->linear_quadtree != NULL)
    {
      printf(" ");
      spamm_ll_print(NULL, node->linear_quadtree);
      printf("%s", empty_header);
    }
    else
    {
      printf("\n");
      printf("%s", empty_header);
    }
    printf("block_dense = %p", (void*) node->block_dense);
    if (node->tier == node->kernel_tier)
    {
      nonzero = 0;

      kernel_block_M = pow(SPAMM_M_CHILD, node->tree_depth-node->kernel_tier)*SPAMM_M_BLOCK;
      kernel_block_N = pow(SPAMM_N_CHILD, node->tree_depth-node->kernel_tier)*SPAMM_N_BLOCK;

      for (i = 0; i < kernel_block_M; ++i) {
        for (j = 0; j < kernel_block_N; ++j)
        {
          if (node->block_dense[spamm_dense_index(i, j, kernel_block_M, kernel_block_N)] != 0.0)
          {
            nonzero++;
          }
        }
      }
      printf(", sparsity = %1.1f%%", (1.0 - (floating_point_t) nonzero / (floating_point_t) (kernel_block_M*kernel_block_N))*100);
      printf(", (1:%u,1:%u) { ", kernel_block_M, kernel_block_N);
      for (i = 0; i < kernel_block_M; ++i) {
        printf(" (%u,1:%u) { ", i+1, kernel_block_N);
        for (j = 0; j < kernel_block_N; ++j)
        {
          printf("%f", node->block_dense[spamm_dense_index(i, j, kernel_block_M, kernel_block_N)]);
          if (j < kernel_block_N-1) { printf(", "); }
        }
        printf(" }");
        if (i < kernel_block_M-1) { printf(", "); }
      }
      printf(" }\n");
    }

    else
    {
      if (node->tier < node->kernel_tier)
      {
        printf(" (above kernel tier)\n");
      }

      else
      {
        printf(" (below kernel tier)\n");
      }
    }
  }
}
