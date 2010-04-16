#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void
spamm_new (const int M, const int N, const int M_block, const int N_block,
    const int M_child, const int N_child, const double threshold,
    struct spamm_t *A)
{
  assert(A != NULL);

  double x, x_M, x_N;

  A->M = M;
  A->N = N;

  A->M_block = M_block;
  A->N_block = N_block;

  A->M_child = M_child;
  A->N_child = N_child;

  /* Pad to powers of M_child x N_child. */
  x_M = log(M >= M_block ? M/(double) M_block : M_block/(double) M)/log(M_child);
  x_N = log(N >= N_block ? N/(double) N_block : N_block/(double) N)/log(N_child);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->M_padded = (int) (M_block*pow(M_child, ceil(x)));
  A->N_padded = (int) (N_block*pow(N_child, ceil(x)));

  //spamm_log("padding matrix from %ix%i to %ix%i\n", __FILE__, __LINE__, M, N, A->M_padded, A->N_padded);

  A->threshold = threshold;

  A->root = NULL;
}
