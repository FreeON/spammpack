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
      return node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_M_BLOCK, SPAMM_N_BLOCK)];
    }

    else if (node->linear_quadtree != NULL)
    {
      /* Read from linear tree. */
      linear_iterator = spamm_ll_iterator_new(node->linear_quadtree);

      LOG2_DEBUG("walking through linear quadtree\n");
      for (linear_tree_node = spamm_ll_iterator_first(linear_iterator); linear_tree_node != NULL; linear_tree_node = spamm_ll_iterator_next(linear_iterator))
      {
        linear_element = linear_tree_node->data;

#ifdef SPAMM_DEBUG
        spamm_int_to_binary(linear_element->index, node->tree_depth*2, binary_string_1);
        LOG_DEBUG("looking at memory starting at %p with index %o (%s)\n", linear_element, linear_element->index, binary_string_1);
#endif

        /* Reconstruct matrix indices from linear index. */
        M_lower = node->M_lower;
        M_upper = node->M_upper;
        N_lower = node->N_lower;
        N_upper = node->N_upper;

        mask_i = 2 << (node->tree_depth-node->tier-1)*2;
        mask_j = 1 << (node->tree_depth-node->tier-1)*2;

#ifdef SPAMM_DEBUG
        spamm_int_to_binary(mask_i, node->tree_depth*2, binary_string_1);
        spamm_int_to_binary(mask_j, node->tree_depth*2, binary_string_2);
        spamm_print_node(node);
        LOG_DEBUG("mask_i = %s, mask_j = %s\n", binary_string_1, binary_string_2);
#endif

        for (l = node->tier; l < node->tree_depth; ++l)
        {
#ifdef SPAMM_DEBUG
          spamm_int_to_binary(linear_element->index & mask_i, node->tree_depth*2, binary_string_1);
          spamm_int_to_binary(linear_element->index & mask_j, node->tree_depth*2, binary_string_2);
          LOG_DEBUG("applying mask_i = %s\n", binary_string_1);
          LOG_DEBUG("applying mask_j = %s\n", binary_string_2);

          LOG_DEBUG("bit_i = %u\n", ((linear_element->index & mask_i) != 0 ? 1 : 0));
          LOG_DEBUG("bit_j = %u\n", ((linear_element->index & mask_j) != 0 ? 1 : 0));
#endif

          M_lower += ((linear_element->index & mask_i) != 0 ? 1 : 0)*(node->M_upper-node->M_lower)/(1 << (l+1-node->tier));
          N_lower += ((linear_element->index & mask_j) != 0 ? 1 : 0)*(node->N_upper-node->N_lower)/(1 << (l+1-node->tier));

#ifdef SPAMM_DEBUG
          LOG_DEBUG("M_lower = %u, N_lower = %u\n", M_lower, N_lower);
#endif

          mask_i = mask_i >> 2;
          mask_j = mask_j >> 2;
        }

        M_upper = M_lower+SPAMM_M_BLOCK;
        N_upper = N_lower+SPAMM_N_BLOCK;

#ifdef SPAMM_DEBUG
        LOG_DEBUG("M_lower = %u, M_upper = %u, N_lower = %u, N_upper = %u\n", M_lower, M_upper, N_lower, N_upper);
#endif

        if (i >= M_lower && i < M_upper && j >= N_lower && j < N_upper)
        {
          spamm_ll_iterator_delete(&linear_iterator);
          return linear_element->block_dense[spamm_dense_index(i-M_lower, j-N_lower, SPAMM_M_BLOCK, SPAMM_N_BLOCK)];
        }
      }
      spamm_ll_iterator_delete(&linear_iterator);
      return 0.0;
    }

    else { return 0.0; }
  }

  else
  {
    /* Recurse down the tree. */
    for (l = 0; l < SPAMM_M_CHILD; ++l) {
      for (k = 0; k < SPAMM_N_CHILD; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/SPAMM_M_CHILD) &&
            i <  (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/SPAMM_M_CHILD) &&
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
    LOG_DEBUG("getting element A(%u,%u)\n", i, j);
    return spamm_get_element(i, j, A->root);
  }

  else { return 0.0; }
}
