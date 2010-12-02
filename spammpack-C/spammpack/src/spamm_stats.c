#include "spamm.h"

/** Return the number of non-zero elements in a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of non-zero elements.
 */
unsigned int
spamm_number_nonzero (const struct spamm_t *A)
{
  unsigned int i, j;
  unsigned int result = 0;

  for (i = 0; i < A->M; i++) {
    for (j = 0; j < A->N; j++)
    {
      if (spamm_get(i, j, A) != 0.0)
      {
        result++;
      }
    }
  }

  return result;
}
