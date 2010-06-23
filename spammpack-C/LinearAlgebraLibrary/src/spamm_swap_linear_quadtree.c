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
  unsigned int i;

  struct spamm_linear_quadtree_t *linear1 = data1;
  struct spamm_linear_quadtree_t *linear2 = data2;

  unsigned int temp_index;
  float_t temp_A;

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

  temp_index = linear1->index;
  linear1->index = linear2->index;
  linear2->index = temp_index;

  for (i = 0; i < linear1->M*linear1->N; ++i)
  {
    temp_A = linear1->block_dense[i];
    linear1->block_dense[i] = linear2->block_dense[i];
    linear2->block_dense[i] = temp_A;
  }
}
