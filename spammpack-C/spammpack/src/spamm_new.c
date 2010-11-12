#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/** Initialize new matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param A The spamm_t matrix.
 */
void
spamm_new (const unsigned int M, const unsigned int N, struct spamm_t *A)
{
  unsigned int depth;
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

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(SPAMM_N_BLOCK) ? log(M) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);
  x_N = (log(N) > log(SPAMM_N_BLOCK) ? log(N) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);

  LOG_DEBUG("x_M = %f, x_N = %f\n", x_M, x_N);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->tree_depth = (unsigned int) ceil(x);

  /* Adjust tree to kernel depth. */
  if (A->tree_depth < SPAMM_KERNEL_DEPTH) { A->tree_depth = SPAMM_KERNEL_DEPTH; }

  /* Set matrix size. */
  A->M = M;
  A->N = N;

  /* Set padded matrix size. */
  A->N_padded = (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, A->tree_depth));

  LOG_DEBUG("creating new SpAMM with M = %u, N = %u\n", M, N);
  LOG_DEBUG("padding to %ux%u\n", A->N_padded, A->N_padded);

  /* Reset the nonzero block counter. */
  A->number_nonzero_blocks = 0;

  /* Set the kernel tier. */
  A->kernel_tier = A->tree_depth-SPAMM_KERNEL_DEPTH;

  A->number_nonzero_blocks = 0;

  A->root = NULL;

  A->norm = 0;

  max_memory = (double) sizeof(struct spamm_t)
    + sizeof(struct spamm_node_t)*pow(SPAMM_N_CHILD*SPAMM_N_CHILD, (double) A->tree_depth)
    + A->N_padded*A->N_padded*sizeof(floating_point_t);

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
