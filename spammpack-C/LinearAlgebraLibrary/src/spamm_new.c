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
spamm_new (const unsigned int M, const unsigned int N, struct spamm_t *A)
{
  double max_memory = 0;

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

  double x, x_M, x_N;

  A->M = M;
  A->N = N;

  /* Pad to powers of M_child x N_child. */
  x_M = fabs(log(M)-log(SPAMM_M_BLOCK))/log(SPAMM_M_CHILD);
  x_N = fabs(log(N)-log(SPAMM_N_BLOCK))/log(SPAMM_N_CHILD);

  LOG_DEBUG("x_M = %f, x_N = %f\n", x_M, x_N);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->number_nonzero_blocks = 0;

  A->M_padded = (int) (SPAMM_M_BLOCK*pow(SPAMM_M_CHILD, ceil(x)));
  A->N_padded = (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, ceil(x)));

  A->tree_depth = (unsigned int) ceil(x);

  /* Set the kernel tier. */
  A->kernel_tier = (A->tree_depth >= SPAMM_KERNEL_DEPTH ? A->tree_depth-SPAMM_KERNEL_DEPTH : 0);

  /* Set the linear tier to depth+1 to indicate that we don't have any linear
   * trees anywhere. */
  A->linear_tier = A->tree_depth+1;

  A->number_nonzero_blocks = 0;

  A->root = NULL;

  max_memory = (double) sizeof(struct spamm_t)
    + sizeof(struct spamm_node_t)*pow(SPAMM_M_CHILD*SPAMM_N_CHILD, (double) A->tree_depth)
    + A->M_padded*A->N_padded*sizeof(floating_point_t);

  if (max_memory < 1024)
  {
    LOG_INFO("max memory for this matrix: %1.2f bytes\n", max_memory);
  }

  else if (max_memory < 1024*1024)
  {
    LOG_INFO("max memory for this matrix: %1.2f kB\n", max_memory/1024.);
  }

  else if (max_memory < 1024*1024*1024)
  {
    LOG_INFO("max memory for this matrix: %1.2f MB\n", max_memory/1024./1024.);
  }

  else
  {
    LOG_INFO("max memory for this matrix: %1.2f GB\n", max_memory/1024./1024./1024.);
  }
}
