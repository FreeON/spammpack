#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

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
    const floating_point_t threshold, struct spamm_t *A)
{
  assert(A != NULL);

  if (M <= 0)
  {
    LOG2_FATAL("M <= 0\n");
    exit(1);
  }

  if (N <= 0)
  {
    LOG2_FATAL("N <= 0\n");
    exit(1);
  }

  if (M_block <= 0)
  {
    LOG2_FATAL("M_block <= 0\n");
    exit(1);
  }

  if (N_block <= 0)
  {
    LOG2_FATAL("N_block <= 0\n");
    exit(1);
  }

  if (M_child <= 0)
  {
    LOG2_FATAL("M_child <= 0\n");
    exit(1);
  }

  if (N_child <= 0)
  {
    LOG2_FATAL("N_child <= 0\n");
    exit(1);
  }

  double x, x_M, x_N;

  A->M = M;
  A->N = N;

  /* Pad to powers of M_child x N_child. */
  x_M = fabs(log(M)-log(M_block))/log(M_child);
  x_N = fabs(log(N)-log(N_block))/log(N_child);

  LOG_DEBUG("x_M = %f, x_N = %f\n", x_M, x_N);

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
