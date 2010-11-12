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
floating_point_t
spamm_get_element (const unsigned int i, const unsigned int j, const struct spamm_node_t *node)
{
  unsigned int l, k;
  unsigned int M_lower, M_upper, N_lower, N_upper;
  unsigned int mask_i, mask_j;

  struct spamm_ll_iterator_t *linear_iterator;
  struct spamm_ll_node_t *linear_tree_node;

  struct spamm_linear_quadtree_t *linear_element;

#ifdef SPAMM_DEBUG
  char binary_string_1[129];
  char binary_string_2[129];
#endif

  assert(node != NULL);

  /* Test whether we are at a leaf node. */
  if (node->tier == node->tree_depth)
  {
    if (node->block_dense != NULL)
    {
      /* Read directly from dense block. */
      return node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_N_BLOCK, SPAMM_N_BLOCK)];
    }

    else { return 0.0; }
  }

  else
  {
    /* Recurse down the tree. */
    for (l = 0; l < SPAMM_N_CHILD; ++l) {
      for (k = 0; k < SPAMM_N_CHILD; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/SPAMM_N_CHILD) &&
            i <  (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/SPAMM_N_CHILD) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/SPAMM_N_CHILD) &&
            j <  (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/SPAMM_N_CHILD))
        {
          if (node->child[l][k] != NULL)
          {
            return spamm_get_element(i, j, node->child[l][k]);
          }

          else { return 0.0; }
        }
      }
    }
  }

  return 0.0;
}

/** Get an element from a matrix tree.
 *
 * @param i Row index.
 * @param j Column index.
 * @param A The matrix tree.
 *
 * @return The matrix element.
 */
floating_point_t
spamm_get (const unsigned int i, const unsigned int j, const struct spamm_t *A)
{
  assert(A != NULL);

  if (i >= A->M) { LOG_FATAL("(i = %i) >= (M = %i)\n", i, A->M); exit(1); }
  if (j >= A->N) { LOG_FATAL("(j = %i) >= (N = %i)\n", j, A->N); exit(1); }

  /* Recurse down to find the element. */
  if (A->root != NULL)
  {
    //LOG_DEBUG("getting element A(%u,%u)\n", i, j);
    return spamm_get_element(i, j, A->root);
  }

  else { return 0.0; }
}
