/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Add two spamm matrices. @f$ A \leftarrow \alpha A + \beta B @f$.
 *
 * @param alpha The factor @f$ \alpha @f$.
 * @param A The matrix A.
 * @param beta The factor @f$ \beta @f$.
 * @param B The matrix B.
 */
void
spamm_hashed_add (const float alpha,
    struct spamm_hashed_t *A,
    const float beta,
    struct spamm_hashed_t *B)
{
  struct spamm_hashtable_t *A_tier_hashtable;
  struct spamm_hashtable_t *B_tier_hashtable;

  assert(A != NULL);
  assert(B != NULL);

  if (A->layout != B->layout)
  {
    printf("[add] inconsisten layout in matrices\n");
    exit(1);
  }

  if (A->M != B->M)
  {
    printf("[add] mismatch of number of rows\n");
    exit(1);
  }

  if (A->N != B->N)
  {
    printf("[add] mismatch of number of columns\n");
    exit(1);
  }

  /* Print out some information. */
  printf("[add] alpha = %e, beta = %e\n", alpha, beta);

  A_tier_hashtable = A->tier_hashtable[0];
  B_tier_hashtable = B->tier_hashtable[0];
}
