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
  unsigned int tier;
  unsigned int *tier_key;

  if (M <= 0)
  {
    fprintf(stderr, "M <= 0\n");
    exit(1);
  }

  if (N <= 0)
  {
    fprintf(stderr, "N <= 0\n");
    exit(1);
  }

  /* Allocate memory. */
  A = (struct spamm_t*) malloc(sizeof(struct spamm_t));

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(SPAMM_N_BLOCK) ? log(M) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);
  x_N = (log(N) > log(SPAMM_N_BLOCK) ? log(N) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);

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

  /* Set the kernel tier. */
  A->kernel_tier = A->depth-SPAMM_KERNEL_DEPTH;

  /* Create the tier hash tables. */
  A->tier_hashtable = g_hash_table_new(g_int_hash, spamm_hash_uint_equal);
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    tier_key = (unsigned int*) malloc(sizeof(unsigned int));
    *tier_key = tier;
    g_hash_table_insert(A->tier_hashtable, tier_key, g_hash_table_new(g_int_hash, spamm_hash_uint_equal));
  }

  return A;
}
