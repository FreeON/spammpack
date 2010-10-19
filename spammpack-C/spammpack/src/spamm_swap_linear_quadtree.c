#include "spamm.h"
#include <stdlib.h>

/** Swap spamm_linear_quadtree_t.
 *
 * @param data1 The first element of type spamm_linear_quadtree_t.
 * @param data2 The second element of type spamm_linear_quadtree_t.
 */
void
spamm_swap_linear_quadtree (void *data1, void *data2)
{
  struct spamm_linear_quadtree_t *linear1 = data1;
  struct spamm_linear_quadtree_t *linear2 = data2;

  if (linear1->M != linear2->M)
  {
    LOG2_FATAL("block dimensions do not match in M\n");
    exit(1);
  }

  if (linear1->N != linear2->N)
  {
    LOG2_FATAL("block dimensions do not match in N\n");
    exit(1);
  }

  spamm_swap_unsigned_int(&linear1->index, &linear2->index);
  spamm_swap_block_dense(linear1->M, linear1->N, linear1->block_dense, linear2->block_dense);
}
