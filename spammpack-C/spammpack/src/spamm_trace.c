/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>

/** The trace of a matrix.
 *
 * @param A The matrix.
 * @param flop The flop count.
 *
 * @return The trace of A.
 */
float
spamm_trace (const struct spamm_matrix_t *const A,
    double *const flop)
{
  unsigned int i[2];
  float trace;

  assert(A != NULL);
  assert(flop != NULL);

  switch(A->number_dimensions)
  {
    case 2:
      if(A->N[0] != A->N[1])
      {
        SPAMM_FATAL("not a square matrix\n");
      }
      for(i[0] = 0, trace = 0; i[0] < A->N[0]; i[0]++)
      {
        i[1] = i[0];
        trace += spamm_get(i, A);
      }
      break;

    default:
      SPAMM_FATAL("not defined\n");
      break;
  }

  *flop += A->N[0];

  return trace;
}
