#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Multiply a matrix node by a scalar. This operation changes the
 * node.
 */
void
spamm_multiply_scalar_node (const floating_point_t alpha, struct spamm_node_t *node)
{
  int i, j;
  struct spamm_ll_iterator_t *iterator;
  struct spamm_ll_node_t *linear_node;
  struct spamm_linear_quadtree_t *linear_element;

  assert(node != NULL);

  if (node->child != NULL)
  {
    for (i = 0; i < node->M_child; ++i) {
      for (j = 0; j < node->N_child; ++j)
      {
        spamm_multiply_scalar_node(alpha, node->child[spamm_dense_index(i, j, node->M_child, node->N_child)]);
      }
    }
  }

  else if (node->block_dense != NULL)
  {
    for (i = 0; i < node->M_block; ++i) {
      for (j = 0; j < node->N_block; ++j)
      {
        node->block_dense[spamm_dense_index(i, j, node->M_block, node->N_block)] *= alpha;
      }
    }
  }

  else if (node->linear_quadtree != NULL)
  {
    iterator = spamm_ll_iterator_new(node->linear_quadtree);
    for (linear_node = spamm_ll_iterator_first(iterator); linear_node != NULL; linear_node = spamm_ll_iterator_next(iterator))
    {
      linear_element = linear_node->data;
      for (i = 0; i < linear_element->M; ++i) {
        for (j = 0; j < linear_element->N; ++j)
        {
          linear_element->block_dense[spamm_dense_index(i, j, linear_element->M, linear_element->N)] *= alpha;
        }
      }
    }
    spamm_ll_iterator_delete(&iterator);
  }
}

/** Multiply the matrix A by a scalar. This operation changes the matrix A.
 *
 * @param alpha The scalar factor to multiply the matrix A by.
 * @param A The matrix.
 */
void
spamm_multiply_scalar (const floating_point_t alpha, struct spamm_t *A)
{
  assert(A != NULL);

  if (A->root != NULL)
  {
    spamm_multiply_scalar_node(alpha, A->root);
  }
}
