/** @file */

#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/** \private Set an element in a matrix.
 *
 * This is the recursive part and is not called directly. Use spamm_set()
 * instead.
 *
 * @param i Row index of element.
 * @param j Column index of element.
 * @param Aij Value to set the matrix element to.
 * @param node The matrix node.
 */
void
spamm_set_element (const int i, const int j, const float_t Aij, struct spamm_node_t *node)
{
  int l, k, m, n;
  struct spamm_node_t *child_node;

  assert(node != NULL);

  /* Check if we are at the block level. */
  if ((node->M_upper-node->M_lower) == node->M_block && (node->N_upper-node->N_lower) == node->N_block)
  {
    if (node->block_dense == NULL)
    {
      spamm_log("bug!\n", __FILE__, __LINE__);
      exit(1);
    }

    node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, node->M_block, node->N_block)] = Aij;
  }

  else
  {
    if (node->child == NULL)
    {
      /* Allocate children nodes. */
      node->child = (struct spamm_node_t**) malloc(sizeof(struct spamm_node_t*)*node->M_child*node->N_child);

      /* Create children nodes. */
      for (l = 0; l < node->M_child; ++l) {
        for (k = 0; k < node->N_child; ++k)
        {
          spamm_new_node(&(node->child[spamm_dense_index(l, k, node->M_child, node->N_child)]));
          child_node = node->child[spamm_dense_index(l, k, node->M_child, node->N_child)];

          child_node->M_lower = node->M_lower+l*(node->M_upper-node->M_lower)/node->M_child;
          child_node->M_upper = node->M_lower+(l+1)*(node->M_upper-node->M_lower)/node->M_child;
          child_node->N_lower = node->N_lower+k*(node->N_upper-node->N_lower)/node->N_child;
          child_node->N_upper = node->N_lower+(k+1)*(node->N_upper-node->N_lower)/node->N_child;

          child_node->M_child = node->M_child;
          child_node->N_child = node->N_child;

          child_node->threshold = node->threshold;

          child_node->M_block = node->M_block;
          child_node->N_block = node->N_block;

          /* Check if we are at the block level. */
          if ((child_node->M_upper-child_node->M_lower) == child_node->M_block && (child_node->N_upper-child_node->N_lower) == child_node->N_block)
          {
            child_node->block_dense = (float_t*) malloc(sizeof(float_t)*child_node->M_block*child_node->N_block);
            for (m = 0; m < child_node->M_block; ++m) {
              for (n = 0; n < child_node->N_block; ++n)
              {
                child_node->block_dense[spamm_dense_index(m, n, child_node->M_block, child_node->N_block)] = 0.0;
              }
            }
          }
        }
      }

      if (node->M_child == 2 && node->N_child == 2)
      {
        /* Z-curve curve ordering. */
        switch (node->ordering)
        {
          case none:
            node->child[spamm_dense_index(0, 0, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(0, 0, node->M_child, node->N_child)]->index = 0;
            node->child[spamm_dense_index(0, 1, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(0, 1, node->M_child, node->N_child)]->index = 1;
            node->child[spamm_dense_index(1, 0, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(1, 0, node->M_child, node->N_child)]->index = 2;
            node->child[spamm_dense_index(1, 1, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(1, 1, node->M_child, node->N_child)]->index = 3;
            break;

          case P:
            node->child[spamm_dense_index(0, 0, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(0, 0, node->M_child, node->N_child)]->index = node->index*4+0;
            node->child[spamm_dense_index(0, 1, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(0, 1, node->M_child, node->N_child)]->index = node->index*4+1;
            node->child[spamm_dense_index(1, 0, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(1, 0, node->M_child, node->N_child)]->index = node->index*4+2;
            node->child[spamm_dense_index(1, 1, node->M_child, node->N_child)]->ordering = P;
            node->child[spamm_dense_index(1, 1, node->M_child, node->N_child)]->index = node->index*4+3;
            break;

          default:
            spamm_log("bug?\n", __FILE__, __LINE__);
            exit(1);
            break;
        }
      }
    }

    /* Recurse. */
    for (l = 0; l < node->M_child; ++l) {
      for (k = 0; k < node->N_child; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/node->M_child) &&
            i < (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/node->M_child) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/node->N_child) &&
            j < (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/node->N_child))
        {
          return spamm_set_element(i, j, Aij, node->child[spamm_dense_index(l, k, node->M_child, node->N_child)]);
        }
      }
    }
  }
}

/** Set an element in a matrix.
 *
 * @param i Row index of element.
 * @param j Column index of element.
 * @param Aij Value to set the matrix element to.
 * @param A The matrix.
 */
int
spamm_set (const int i, const int j, const float_t Aij, struct spamm_t *A)
{
  int l, k;
  int result = 0;

  assert(A != NULL);

  if (i < 0)     { spamm_log("(i = %i) < 0\n",  __FILE__, __LINE__, i); exit(1); }
  if (j < 0)     { spamm_log("(j = %i) < 0\n",  __FILE__, __LINE__, j); exit(1); }
  if (i >= A->M) { spamm_log("(i = %i) >= (M = %i)\n", __FILE__, __LINE__, i, A->M); exit(1); }
  if (j >= A->N) { spamm_log("(j = %i) >= (N = %i)\n", __FILE__, __LINE__, j, A->N); exit(1); }

  if (fabs(Aij) > A->threshold)
  {
    /* Recursively find the leaf node that stores this element. */
    if (A->root == NULL)
    {
      spamm_new_node(&(A->root));
      A->root->M_lower = 0;
      A->root->M_upper = A->M_padded;
      A->root->N_lower = 0;
      A->root->N_upper = A->N_padded;

      A->root->M_child = A->M_child;
      A->root->N_child = A->N_child;

      A->root->threshold = A->threshold;

      A->root->M_block = A->M_block;
      A->root->N_block = A->N_block;

      /* Check if we are at the block level. */
      if ((A->root->M_upper-A->root->M_lower) == A->root->M_block && (A->root->N_upper-A->root->N_lower) == A->root->N_block)
      {
        A->root->block_dense = (float_t*) malloc(sizeof(float_t)*A->root->M_block*A->root->N_block);
        for (l = 0; l < A->root->M_block; ++l) {
          for (k = 0; k < A->root->N_block; ++k)
          {
            A->root->block_dense[spamm_dense_index(l, k, A->root->M_block, A->root->N_block)] = 0.0;
          }
        }
      }
    }

    spamm_set_element(i, j, Aij, A->root);
  }

  else
  {
    result = -1;
  }

  return result;
}
