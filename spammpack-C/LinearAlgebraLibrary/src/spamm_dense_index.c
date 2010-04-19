#include "spamm.h"
#include <assert.h>

int
spamm_dense_index (const int i, const int j, const int M, const int N)
{
  /* Row-Major format. */
  //return i*N+j;

  /* Column-Major format.
   *
   * When using dgemm, we need column-major storage.
   */
  return i+j*M;
}
