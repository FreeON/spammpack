#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

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

  /* We need at least one level. */
  //if (ceil(x) < 1.0) { x = 1.0; }

  A->M_padded = (int) (M_block*pow(M_child, ceil(x)));
  A->N_padded = (int) (N_block*pow(N_child, ceil(x)));

  //spamm_log("M = %i, N = %i\n", __FILE__, __LINE__, M, N);
  //spamm_log("M_block = %i, N_block = %i\n", __FILE__, __LINE__, M_block, N_block);
  //spamm_log("M_child = %i, N_child = %i\n", __FILE__, __LINE__, M_child, N_child);
  //spamm_log("x_M = %e, x_N = %e\n", __FILE__, __LINE__, x_M, x_M);
  //spamm_log("padding matrix from %ix%i to %ix%i\n", __FILE__, __LINE__, M, N, A->M_padded, A->N_padded);

  A->threshold = threshold;

  A->root = NULL;
}
