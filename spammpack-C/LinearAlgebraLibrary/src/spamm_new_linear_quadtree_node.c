#include "spamm.h"

/** Create new linear quadtree block.
 *
 * @param M The number of rows of the dense matrix block.
 * @param N The number of columns of the dense matrix block.
 * @param memory The memory to allocate the linear tree node in.
 *
 * @return A new node in the linear quadtree.
 */
struct spamm_linear_quadtree_t*
spamm_new_linear_quadtree_node (const unsigned int M, const unsigned int N,
    struct spamm_mm_t *memory)
{
  struct spamm_linear_quadtree_t *node;

  node = (struct spamm_linear_quadtree_t*) spamm_mm_allocate(sizeof(struct spamm_linear_quadtree_t)+M*N*sizeof(float_t)+8, memory);
  node->block_dense = (float_t*) (((void*) node)+sizeof(struct spamm_linear_quadtree_t));
  node->M = M;
  node->N = N;

  return node;
}
