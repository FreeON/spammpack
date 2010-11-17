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
 *
 * @return The square of the Frobenius norm of this node.
 */
floating_point_t
spamm_set_element (const unsigned int i, const unsigned int j, const floating_point_t Aij, struct spamm_node_t *node)
{
  int l, k, m, n;
  unsigned int kernel_block_N;
  struct spamm_node_t *child_node;
  char binary_string[1000];

  assert(node != NULL);

  /* Check if we are at the block level. */
  if (node->tier == node->tree_depth)
  {
    if (node->block_dense == NULL)
    {
      LOG2_FATAL("bug!\n");
      exit(1);
    }

    /* Set matrix element. */
    node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_N_BLOCK, SPAMM_N_BLOCK)] = Aij;

    /* Set the dilated matrix element. */
    node->block_dense_dilated[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_N_BLOCK, SPAMM_N_BLOCK)*4+0] = Aij;
    node->block_dense_dilated[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_N_BLOCK, SPAMM_N_BLOCK)*4+1] = Aij;
    node->block_dense_dilated[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_N_BLOCK, SPAMM_N_BLOCK)*4+2] = Aij;
    node->block_dense_dilated[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_N_BLOCK, SPAMM_N_BLOCK)*4+3] = Aij;

    /* Recalculate the norm. */
    node->norm2 += Aij*Aij;
    node->norm = sqrt(node->norm2);

    return node->norm2;
  }

  else
  {
    /* Recurse. */
    for (l = 0; l < SPAMM_N_CHILD; l++) {
      for (k = 0; k < SPAMM_N_CHILD; k++)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/SPAMM_N_CHILD) &&
            i <  (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/SPAMM_N_CHILD) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/SPAMM_N_CHILD) &&
            j <  (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/SPAMM_N_CHILD))
        {
          if (node->child[l][k] == NULL)
          {
            /* Create new child node. */
            node->child[l][k] = spamm_new_childnode(node->tier+1, node->tree_depth,
                node->M_lower+l*(node->M_upper-node->M_lower)/SPAMM_N_CHILD,
                node->M_lower+(l+1)*(node->M_upper-node->M_lower)/SPAMM_N_CHILD,
                node->N_lower+k*(node->N_upper-node->N_lower)/SPAMM_N_CHILD,
                node->N_lower+(k+1)*(node->N_upper-node->N_lower)/SPAMM_N_CHILD,
                node->M_lower_kernel_tier, node->M_upper_kernel_tier,
                node->N_lower_kernel_tier, node->N_upper_kernel_tier,
                node->kernel_tier, node->block_dense, node->block_dense_dilated);

            if (SPAMM_N_CHILD == 2)
            {
              /* Set the linear index on child node. */
              node->child[l][k]->index = (node->index << 2) | (l << 1) | k;
              spamm_int_to_binary(node->child[l][k]->index, 2*(node->tier+1), binary_string);
              LOG_DEBUG("setting linear index of child[%u][%u] to %s\n", l, k, binary_string);
            }
          }

          node->norm2 -= node->child[l][k]->norm2;
          node->norm2 += spamm_set_element(i, j, Aij, node->child[l][k]);
          node->norm = sqrt(node->norm2);
        }
      }
    }

    return node->norm2;
  }
}

/** Set an element in a matrix.
 *
 * @param i Row index of element.
 * @param j Column index of element.
 * @param Aij Value to set the matrix element to.
 * @param A The matrix.
 *
 * @return #SPAMM_RESULT_OK indicates that everything went fine.
 * @return #SPAMM_RESULT_ZERO_ELEMENT indicates that the matrix element was
 *         zero and therefore not stored.
 */
int
spamm_set (const unsigned int i, const unsigned int j, const floating_point_t Aij, struct spamm_t *A)
{
  int l, k;
  int result = SPAMM_RESULT_OK;

  assert(A != NULL);

  if (i >= A->M) { LOG_FATAL("(i = %i) >= (M = %i)\n", i, A->M); exit(1); }
  if (j >= A->N) { LOG_FATAL("(j = %i) >= (N = %i)\n", j, A->N); exit(1); }

  /* If the value is zero, we don't have to store it. */
  if (Aij != 0.0)
  {
    /* Recursively find the leaf node that stores this element. */
    if (A->root == NULL)
    {
      A->root = spamm_new_childnode(0, A->tree_depth,
          0, A->N_padded, 0, A->N_padded,
          0, 0, 0, 0, A->kernel_tier, NULL, NULL);
    }

    spamm_set_element(i, j, Aij, A->root);

    /* Set norm. */
    A->norm = A->root->norm;
  }

  else
  {
    result = SPAMM_RESULT_ZERO_ELEMENT;
  }

  return result;
}
