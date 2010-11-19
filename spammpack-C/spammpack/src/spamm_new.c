#include "spamm.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Initialize new matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 * @param A The spamm_t matrix.
 */
struct spamm_t *
spamm_new (const unsigned int M, const unsigned int N)
{
  struct spamm_t *A;
  double x, x_M, x_N;

  if (M <= 0)
  {
    printf("M <= 0\n");
    exit(1);
  }

  if (N <= 0)
  {
    printf("N <= 0\n");
    exit(1);
  }

  /* Allocate memory. */
  A = (struct spamm_t*) malloc(sizeof(struct spamm_t));

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(SPAMM_N_BLOCK) ? log(M) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);
  x_N = (log(N) > log(SPAMM_N_BLOCK) ? log(N) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);

  printf("x_M = %f, x_N = %f\n", x_M, x_N);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  A->depth = (unsigned int) ceil(x);

  /* Adjust tree to kernel depth. */
  if (A->depth < SPAMM_KERNEL_DEPTH) { A->depth = SPAMM_KERNEL_DEPTH; }

  /* Set matrix size. */
  A->M = M;
  A->N = N;

  /* Set padded matrix size. */
  A->N_padded = (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, A->depth));

  printf("created %ux%u matrix, padded to %ux%u\n", M, N, A->N_padded, A->N_padded);

  /* Set the kernel tier. */
  A->kernel_tier = A->depth-SPAMM_KERNEL_DEPTH;

  /* Create the tier hash tables. */
  A->tier = g_hash_table_new(g_int_hash, g_int_equal);

  return A;
}
