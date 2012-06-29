/** @file */

#include "spamm.h"

#include <assert.h>

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
  assert(A != NULL);
  assert(B != NULL);


}
