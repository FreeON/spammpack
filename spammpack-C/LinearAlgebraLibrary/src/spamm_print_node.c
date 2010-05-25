/** @file */

#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** \private Convert an integer to string in binary.
 *
 * @param integer The integer to convert.
 * @param width The width of the string, i.e. how many bits should be printed.
 * @param binary_string The string.
 */
void
spamm_int_to_binary (const unsigned int integer, const int width, char *binary_string)
{
  int i;
  unsigned int mask = 1;

  if (width == 0)
  {
    binary_string[0] = '0';
    binary_string[1] = '\0';
  }

  else
  {
    for (i = 0; i < width; ++i)
    {
      binary_string[i] = '0';
    }
    binary_string[width] = '\0';

    for (i = 0; i < width; ++i)
    {
      if (mask & integer) { binary_string[width-1-i] = '1'; }
      mask = mask << 1;
    }
  }
}

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
    printf("linear = %i, ", node->linear_tiers);
    printf("M_lower = %i, M_upper = %i, ", node->M_lower, node->M_upper);
    printf("N_lower = %i, N_upper = %i, ", node->N_lower, node->N_upper);
    printf("M_block = %i, N_block = %i, ", node->M_block, node->N_block);
    printf("M_child = %i, N_child = %i\n", node->M_child, node->N_child);
    printf("%s", empty_header);
    printf("threshold = %7.1e, ", node->threshold);
    binary_string = (char*) malloc(sizeof(char)*(node->tier*2+1));
    spamm_int_to_binary(node->index, node->tier*2, binary_string);
    printf("index = 0b%s (%u), ", binary_string, node->index);
    free(binary_string);
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
      printf(" }\n");
      printf("%s", empty_header);
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
      printf(", sparsity = %1.1f%%", (1.0 - (float_t) nonzero / (float_t) (node->M_block*node->N_block))*100);
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
