/** @file */

#include "spamm.h"
#include "spamm_types_private.h"

float
spamm_trace (const struct spamm_matrix_t *const A)
{
  unsigned int i[2];
  float trace;

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

  SPAMM_WARN("trace = %e\n", trace);
}
