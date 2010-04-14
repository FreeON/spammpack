#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void
spamm_new (const int M, const int N, const int M_block, const int N_block, struct spamm_t *A)
{
  assert(A != NULL);
  assert(M > 0 && N > 0);

  double x, x_M, x_N;

  A->M = M;
  A->N = N;

  A->M_block = M_block;
  A->N_block = N_block;

  /* Pad to even powers of 2. */
  x_M = log(M/(double) M_block)/log(2);
  x_N = log(N/(double) N_block)/log(2);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->M_padded = (int) (M_block*pow(2, ceil(x)));
  A->N_padded = (int) (N_block*pow(2, ceil(x)));

  spamm_log("padding matrix from %ix%i to %ix%i\n", __FILE__, __LINE__, M, N, A->M_padded, A->N_padded);

  A->root = NULL;
}
