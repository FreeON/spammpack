#include "spamm.h"
#include <stdio.h>
#include <stdlib.h>

void
spamm_print_node (const struct spamm_node_t *node)
{
  printf("node:\n");
  printf("  block dimensions: %ix%i\n", node->M_block, node->N_block);
  printf("  node coverage: [%i-%i[ x [%i-%i[\n", node->M_lower, node->M_upper, node->N_lower, node->N_upper);
}
