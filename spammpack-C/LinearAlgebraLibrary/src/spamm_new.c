#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/** Initialize new node of matrix.
 *
 * @param node The spamm_node_t node to initialize.
 */
void
spamm_new_node (struct spamm_node_t **node)
{
  *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  (*node)->tier = 0;
  (*node)->tree_depth = 0;
  (*node)->linear_tier = 0;

  (*node)->M_upper = 0;
  (*node)->M_lower = 0;
  (*node)->N_upper = 0;
  (*node)->N_lower = 0;

  (*node)->M_child = 0;
  (*node)->N_child = 0;

  (*node)->M_block = 0;
  (*node)->N_block = 0;

  (*node)->threshold = 0.0;

  (*node)->index = 0;

  (*node)->previous_i = NULL;
  (*node)->next_i = NULL;
  (*node)->previous_j = NULL;
  (*node)->next_j = NULL;

  (*node)->ordering = none;

  (*node)->block_loaded_in_GPU = 0;

  (*node)->child = NULL;

  (*node)->linear_quadtree = NULL;
  (*node)->block_dense = NULL;
}

/** Initialize new matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param M_block Number of rows of matrix block in spamm_t.
 * @param N_block Number of columns of matrix block in spamm_t.
 * @param M_child Number of rows of children array in spamm_node_t.
 * @param N_child Number of columns of children array in spamm_node_t.
 * @param threshold Threshold below which matrix elements are considered zero.
 * @param A The spamm_t matrix.
 */
void
spamm_new (const unsigned int M, const unsigned int N,
    const unsigned int M_block, const unsigned int N_block,
    const unsigned int M_child, const unsigned int N_child,
    const float_t threshold, struct spamm_t *A)
{
  assert(A != NULL);

  if (M <= 0)
  {
    spamm_log("M <= 0\n", __FILE__, __LINE__);
    exit(1);
  }

  if (N <= 0)
  {
    spamm_log("N <= 0\n", __FILE__, __LINE__);
    exit(1);
  }

  if (M_block <= 0)
  {
    spamm_log("M_block <= 0\n", __FILE__, __LINE__);
    exit(1);
  }

  if (N_block <= 0)
  {
    spamm_log("N_block <= 0\n", __FILE__, __LINE__);
    exit(1);
  }

  if (M_child <= 0)
  {
    spamm_log("M_child <= 0\n", __FILE__, __LINE__);
    exit(1);
  }

  if (N_child <= 0)
  {
    spamm_log("N_child <= 0\n", __FILE__, __LINE__);
    exit(1);
  }

  double x, x_M, x_N;

  A->M = M;
  A->N = N;

  /* Pad to powers of M_child x N_child. */
  x_M = fabs(log(M)-log(M_block))/log(M_child);
  x_N = fabs(log(N)-log(N_block))/log(N_child);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->number_nonzero_blocks = 0;

  A->M_padded = (int) (M_block*pow(M_child, ceil(x)));
  A->N_padded = (int) (N_block*pow(N_child, ceil(x)));

  A->M_block = M_block;
  A->N_block = N_block;

  A->M_child = M_child;
  A->N_child = N_child;

  A->threshold = threshold;

  A->tree_depth = (unsigned int) ceil(x);

  /* Set the linear tier to depth+1 to indicate that we don't have any linear
   * trees anywhere. */
  A->linear_tier = A->tree_depth+1;

  A->number_nonzero_blocks = 0;

  A->root = NULL;
}
