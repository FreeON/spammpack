#include "spamm.h"
#include <assert.h>
#include <stdlib.h>

/** \private Get an element from a matrix tree.
 *
 * This is the recursive part and is not called directly. Use spamm_get().
 *
 * @param i Row index.
 * @param j Column index.
 * @param node The matrix node.
 *
 * @return The matrix element.
 */
float_t
spamm_get_element (const unsigned int i, const unsigned int j, const struct spamm_node_t *node)
{
  unsigned int l, k;
  struct spamm_ll_node_t *linear_tree_node_outer, *linear_tree_node_inner;
  struct spamm_linear_quadtree_t *linear_element;
  float_t result = 0.0;

  assert(node != NULL);

  /* Test whether we are at a leaf node. */
  if (node->child == NULL)
  {
    if (node->block_dense != NULL)
    {
      /* Read directly from dense block. */
      return node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, node->M_block, node->N_block)];
    }

    else if (node->linear_quadtree != NULL)
    {
      /* Read from linear tree. */
      for (linear_tree_node_outer = node->linear_quadtree->first; linear_tree_node_outer != NULL; linear_tree_node_outer = linear_tree_node_outer->next) {
        for (linear_tree_node_inner = ((struct spamm_mm_t*) linear_tree_node_outer->data)->allocated_start.first; linear_tree_node_inner != NULL; linear_tree_node_inner = linear_tree_node_inner->next)
        {
          linear_element = linear_tree_node_inner->data;
          LOG("looking at memory starting at %p with index %o\n", linear_element, linear_element->index);

          /* [FIXME] We need to translate (i,j) into a linear index by bit
           * interleaving. Then we can identify the linear_element we really
           * want and what matrix element inside its dense matrix we need.
           */
          for (l = node->tier; l < node->tree_depth; ++l)
          {
          }
        }
      }
    }

    else { return 0.0; }
  }

  else
  {
    /* Recurse down the tree. */
    for (l = 0; l < node->M_child; ++l) {
      for (k = 0; k < node->N_child; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/node->M_child) &&
            i <  (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/node->M_child) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/node->N_child) &&
            j <  (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/node->N_child))
        {
          result = spamm_get_element(i, j, node->child[spamm_dense_index(l, k, node->M_child, node->N_child)]);
        }
      }
    }
  }

  return result;
}

/** Get an element from a matrix tree.
 *
 * @param i Row index.
 * @param j Column index.
 * @param A The matrix tree.
 *
 * @return The matrix element.
 */
float_t
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A)
{
  assert(A != NULL);

  if (i >= A->M) { LOG("(i = %i) >= (M = %i)\n", i, A->M); exit(1); }
  if (j >= A->N) { LOG("(j = %i) >= (N = %i)\n", j, A->N); exit(1); }

  /* Recurse down to find the element. */
  if (A->root != NULL)
  {
    return spamm_get_element(i, j, A->root);
  }

  else { return 0.0; }
}
