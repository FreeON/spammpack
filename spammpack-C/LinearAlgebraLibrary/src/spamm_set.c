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
spamm_set_element (const unsigned int i, const unsigned int j, const floating_point_t Aij, struct spamm_node_t *node)
{
  int l, k, m, n;
  struct spamm_node_t *child_node;
  unsigned int kernel_block_M, kernel_block_N;

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
    node->block_dense[spamm_dense_index(i-node->M_lower, j-node->N_lower, SPAMM_M_BLOCK, SPAMM_N_BLOCK)] = Aij;
  }

  else
  {
    /* Recurse. */
    for (l = 0; l < SPAMM_M_CHILD; ++l) {
      for (k = 0; k < SPAMM_N_CHILD; ++k)
      {
        if (i >= (node->M_lower+(node->M_upper-node->M_lower)*l/SPAMM_M_CHILD) &&
            i < (node->M_lower+(node->M_upper-node->M_lower)*(l+1)/SPAMM_M_CHILD) &&
            j >= (node->N_lower+(node->N_upper-node->N_lower)*k/SPAMM_N_CHILD) &&
            j < (node->N_lower+(node->N_upper-node->N_lower)*(k+1)/SPAMM_N_CHILD))
        {
          if (node->child[l][k] == NULL)
          {
            /* Create new child node. */
            node->child[l][k] = spamm_new_node();
            child_node = node->child[l][k];

            child_node->tier = node->tier+1;
            child_node->tree_depth = node->tree_depth;

            child_node->M_lower = node->M_lower+l*(node->M_upper-node->M_lower)/SPAMM_M_CHILD;
            child_node->M_upper = node->M_lower+(l+1)*(node->M_upper-node->M_lower)/SPAMM_M_CHILD;
            child_node->N_lower = node->N_lower+k*(node->N_upper-node->N_lower)/SPAMM_N_CHILD;
            child_node->N_upper = node->N_lower+(k+1)*(node->N_upper-node->N_lower)/SPAMM_N_CHILD;

            child_node->linear_tier = node->linear_tier;
            child_node->kernel_tier = node->kernel_tier;

            /* Check if we are at the kernel level. */
            if (child_node->tier == child_node->kernel_tier)
            {
              /* Allocate contiguous matrix block. */
              kernel_block_M = pow(SPAMM_M_CHILD, child_node->tree_depth-child_node->kernel_tier)*SPAMM_M_BLOCK;
              kernel_block_N = pow(SPAMM_N_CHILD, child_node->tree_depth-child_node->kernel_tier)*SPAMM_N_BLOCK;

              child_node->block_dense = (floating_point_t*) spamm_allocate(sizeof(floating_point_t)*kernel_block_M*kernel_block_N);
              for (m = 0; m < kernel_block_M; ++m) {
                for (n = 0; n < kernel_block_N; ++n)
                {
                  child_node->block_dense[spamm_dense_index(m, n, kernel_block_M, kernel_block_N)] = 0.0;
                }
              }
            }

            else if (child_node->tier > child_node->kernel_tier)
            {
              /* Point into the contiguous matrix block. */
              kernel_block_M = pow(SPAMM_M_CHILD, child_node->tree_depth-child_node->tier)*SPAMM_M_BLOCK;
              kernel_block_N = pow(SPAMM_N_CHILD, child_node->tree_depth-child_node->tier)*SPAMM_N_BLOCK;

              child_node->block_dense = node->block_dense+kernel_block_M*kernel_block_N*(SPAMM_N_CHILD*l+k);
            }
          }

          return spamm_set_element(i, j, Aij, node->child[l][k]);
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
 *
 * @return #SPAMM_RESULT_OK indicates that everything went fine.
 * @return #SPAMM_RESULT_ZERO_ELEMENT indicates that the matrix element was
 *         zero and therefore not stored.
 */
int
spamm_set (const unsigned int i, const unsigned int j, const floating_point_t Aij, struct spamm_t *A)
{
  int l, k;
  unsigned int kernel_block_M, kernel_block_N;
  int result = SPAMM_RESULT_OK;

  assert(A != NULL);

  if (i >= A->M) { LOG_FATAL("(i = %i) >= (M = %i)\n", i, A->M); exit(1); }
  if (j >= A->N) { LOG_FATAL("(j = %i) >= (N = %i)\n", j, A->N); exit(1); }

  /* If the value is zero, we don't have to store it. */
  if (Aij != 0.0)
  {
    //printf("setting A(%u,%u) to %f\n", i+1, j+1, Aij);

    /* Recursively find the leaf node that stores this element. */
    if (A->root == NULL)
    {
      A->root = spamm_new_node();

      A->root->tree_depth = A->tree_depth;

      A->root->M_lower = 0;
      A->root->M_upper = A->M_padded;
      A->root->N_lower = 0;
      A->root->N_upper = A->N_padded;

      A->root->linear_tier = A->linear_tier;
      A->root->kernel_tier = A->kernel_tier;

      /* Check if we are at the kernel level. */
      if (A->root->tier == A->root->kernel_tier)
      {
        kernel_block_M = pow(SPAMM_M_CHILD, A->root->tree_depth-A->root->kernel_tier)*SPAMM_M_BLOCK;
        kernel_block_N = pow(SPAMM_N_CHILD, A->root->tree_depth-A->root->kernel_tier)*SPAMM_N_BLOCK;

        /* Reset newly allocated block to zero. */
        A->root->block_dense = (floating_point_t*) spamm_allocate(sizeof(floating_point_t)*kernel_block_M*kernel_block_N);
        for (l = 0; l < kernel_block_M; ++l) {
          for (k = 0; k < kernel_block_N; ++k)
          {
            A->root->block_dense[spamm_dense_index(l, k, kernel_block_M, kernel_block_N)] = 0.0;
          }
        }
      }
    }

    spamm_set_element(i, j, Aij, A->root);
    //spamm_print_tree(A);
  }

  else
  {
    result = SPAMM_RESULT_ZERO_ELEMENT;
  }

  return result;
}
