#include "spamm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/** Initialize new matrix object.
 *
 * @param M Number of rows of dense input matrix.
 * @param N Number of columns of dense input matrix.
 *
 * @return A pointer to the matrix.
 */
struct spamm_t *
spamm_new (const unsigned int M, const unsigned int N, const enum spamm_layout_t layout)
{
  struct spamm_t *A;
  double x, x_M, x_N;
  unsigned int tier;

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

  /* Set the layout. */
  switch (layout)
  {
    case row_major:
    case column_major:
    case Z_curve:
      A->layout = layout;
      break;

    default:
      fprintf(stderr, "[spamm new] unknown layout (%i)\n", layout);
      exit(1);
      break;
  }

  /* Pad to powers of M_child x N_child. */
  x_M = (log(M) > log(SPAMM_N_BLOCK) ? log(M) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);
  x_N = (log(N) > log(SPAMM_N_BLOCK) ? log(N) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  /* The ceil() function can lead to a depth that is one tier too large
   * because of numerical errors in the calculation of x. We need to check
   * whether the depth is appropriate.
   */
  A->depth = (unsigned int) ceil(x);

  /* Double check depth. */
  if (A->depth >= 1 && ((int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, A->depth-1)) >= M && (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, A->depth-1)) >= N))
  {
    (A->depth)--;
  }

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
  A->tier_hashtable = (struct spamm_hashtable_t**) malloc(sizeof(struct spamm_hashtable_t*)*(A->kernel_tier+1));
  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    A->tier_hashtable[tier] = spamm_hashtable_new();
  }

  return A;
}
