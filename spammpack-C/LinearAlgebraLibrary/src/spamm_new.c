/** @file */

#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/** Initialize new matrix object.
 */
void
spamm_new (const int M, const int N, const int M_block, const int N_block,
    const int M_child, const int N_child, const float_t threshold,
    struct spamm_t *A)
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

  A->M_block = M_block;
  A->N_block = N_block;

  A->M_child = M_child;
  A->N_child = N_child;

  /* Pad to powers of M_child x N_child. */
  x_M = fabs(log(M)-log(M_block))/log(M_child);
  x_N = fabs(log(N)-log(N_block))/log(N_child);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->tree_depth = (unsigned int) ceil(x);

  A->linear_tiers = 0;
  A->number_nonzero_blocks = 0;

  A->M_padded = (int) (M_block*pow(M_child, ceil(x)));
  A->N_padded = (int) (N_block*pow(N_child, ceil(x)));

  A->threshold = threshold;

  A->number_nonzero_blocks = 0;

  A->root = NULL;
}

/** Initialize new node of matrix.
 */
void
spamm_new_node (struct spamm_node_t **node)
{
  *node = (struct spamm_node_t*) malloc(sizeof(struct spamm_node_t));

  (*node)->tier = 0;
  (*node)->linear_tiers = 0;

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
  (*node)->ordering = none;

  (*node)->block_loaded_in_GPU = 0;

  (*node)->child = NULL;
  (*node)->block_dense = NULL;
}
