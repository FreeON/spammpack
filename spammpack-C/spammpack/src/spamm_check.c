#include "spamm.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

/** Check the internal consistency of a matrix.
 *
 * @param A The matrix to check
 *
 * @return The following error codes are returned:
 *   - SPAMM_OK - The matrix is consistent.
 *   - SPAMM_ERROR - Something is not consistent.
 */
int
spamm_check (const struct spamm_t *A)
{
  int result = SPAMM_OK;
  unsigned int depth;
  unsigned int N_padded;
  unsigned int tier;
  unsigned int reverse_tier;
  float x_M, x_N, x;
  GHashTable *hashtable;

  assert(A != NULL);

  /* Calculate padding and depth of matrix based on values stored in M and N.
   */
  x_M = (log(A->M) > log(SPAMM_N_BLOCK) ? log(A->M) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);
  x_N = (log(A->N) > log(SPAMM_N_BLOCK) ? log(A->N) - log(SPAMM_N_BLOCK) : 0)/log(SPAMM_N_CHILD);

  if (x_M > x_N) { x = x_M; }
  else           { x = x_N; }

  depth = (unsigned int) ceil(x);

  if (A->depth != depth)
  {
    printf("depth incorrect, should be %u, but is %u\n", depth, A->depth);
    return SPAMM_ERROR;
  }

  N_padded = (int) (SPAMM_N_BLOCK*pow(SPAMM_N_CHILD, depth));

  if (A->N_padded != N_padded)
  {
    printf("padded matrix dimensions incorrect, should be %u, but is %u\n", N_padded, A->N_padded);
    return SPAMM_ERROR;
  }

  if (A->kernel_tier != depth-SPAMM_KERNEL_DEPTH)
  {
    printf("kernel tier incorrect, should be %u, but is %u\n", depth-SPAMM_KERNEL_DEPTH, A->kernel_tier);
    return SPAMM_ERROR;
  }

  /* Check whether there are tier hashtables for every tier. */
  if (A->tier_hashtable == NULL)
  {
    printf("no tier hashtable\n");
    return SPAMM_ERROR;
  }

  for (tier = 0; tier <= A->kernel_tier; tier++)
  {
    if ((hashtable = g_hash_table_lookup(A->tier_hashtable, &tier)) == NULL)
    {
      printf("missing tier hashtable for tier %u\n", tier);
      return SPAMM_ERROR;
    }
  }

  /* Check each node. */
  for (tier = A->kernel_tier+1; tier >= 1; tier--)
  {
    /* Get tier hashtable. */
    reverse_tier = tier-1;
    hashtable = g_hash_table_lookup(A->tier_hashtable, &reverse_tier);

    /* Verify consistency of 2D and 3D linear indices. */
  }

  return result;
}
