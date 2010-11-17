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
  char binary_string[2000];
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
    printf("A(%u:%u,%u:%u), ", node->M_lower+1, node->M_upper, node->N_lower+1, node->N_upper);
    printf("norm = %f, ", node->norm);
    printf("M_lower = %i, M_upper = %i, ", node->M_lower, node->M_upper);
    printf("N_lower = %i, N_upper = %i, ", node->N_lower, node->N_upper);
    printf("\n");
    printf("%s", empty_header);
    spamm_int_to_binary(node->index_2D, node->tier*2, binary_string);
    printf("index_2D = 0b%s (%u), ", binary_string, node->index_2D);
    spamm_int_to_binary(node->index_3D_column, node->tier*3, binary_string);
    printf("index_3D_column = 0b%s (%u), ", binary_string, node->index_3D_column);
    spamm_int_to_binary(node->index_3D_row, node->tier*3, binary_string);
    printf("index_3D_row = 0b%s (%u), ", binary_string, node->index_3D_row);
    printf("child = %p, ", (void*) node->child);
    if (node->child != NULL)
    {
      printf("(1:%u,1:%u) { ", SPAMM_N_CHILD, SPAMM_N_CHILD);
      for (i = 0; i < SPAMM_N_CHILD; ++i) {
        printf("(%u,1:%u) { ", i+1, SPAMM_N_CHILD);
        for (j = 0; j < SPAMM_N_CHILD; ++j)
        {
          printf("%p", (void*) node->child[i][j]);
          if (j < SPAMM_N_CHILD-1) { printf(", "); }
        }
        printf(" }");
        if (i < SPAMM_N_CHILD-1) { printf(", "); }
      }
      printf(" }\n");
      printf("%s", empty_header);
    }
    printf("block_dense = %p", (void*) node->block_dense);
    if (node->tier == node->kernel_tier)
    {
      nonzero = 0;
      for (i = 0; i < SPAMM_N_KERNEL; ++i) {
        for (j = 0; j < SPAMM_N_KERNEL; ++j)
        {
          if (node->block_dense[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)] != 0.0)
          {
            nonzero++;
          }
        }
      }
      printf(", sparsity = %1.1f%%", (1.0 - (floating_point_t) nonzero / (floating_point_t) (SPAMM_N_KERNEL*SPAMM_N_KERNEL))*100);
      printf(", (1:%u,1:%u) { ", SPAMM_N_KERNEL, SPAMM_N_KERNEL);
      for (i = 0; i < SPAMM_N_KERNEL; ++i) {
        printf(" (%u,1:%u) { ", i+1, SPAMM_N_KERNEL);
        for (j = 0; j < SPAMM_N_KERNEL; ++j)
        {
          printf("%0.3f", node->block_dense[spamm_dense_index(i, j, SPAMM_N_KERNEL, SPAMM_N_KERNEL)]);
          if (j < SPAMM_N_KERNEL-1) { printf(", "); }
        }
        printf(" }");
        if (i < SPAMM_N_KERNEL-1) { printf(", "); }
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
