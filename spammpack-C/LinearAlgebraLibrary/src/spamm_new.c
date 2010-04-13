#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void
spamm_new (const int M, const int N, const int M_block, const int N_block, struct spamm_t *A)
{
  assert(A != NULL);
  assert(M > 0 && N > 0);

  A->M = M;
  A->N = N;

  A->M_block = M_block;
  A->N_block = N_block;

  /* Pad to even powers of 2. */
  A->M_padded = (int) pow(2, ceil(log(M/M_block)/log(2)));
  A->N_padded = (int) pow(2, ceil(log(N/N_block)/log(2)));

  spamm_log("padding matrix from %ix%i to %ix%i\n", __FILE__, __LINE__, M, N, A->M_padded, A->N_padded);

  A->root = NULL;
}
