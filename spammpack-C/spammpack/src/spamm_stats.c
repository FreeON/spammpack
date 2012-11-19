#include "spamm.h"
#include "spamm_types_private.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

/** Return the number of non-zero elements in a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of non-zero elements.
 */
unsigned int
spamm_number_nonzero (const struct spamm_matrix_t *A)
{
  unsigned int *i;
  unsigned int result = 0;

  i = calloc(A->number_dimensions, sizeof(unsigned int));

  switch (A->number_dimensions)
  {
    case 1:
      for (i[0] = 0; i[0] < A->N[0]; i[0]++)
      {
        if (spamm_get(i, A) != 0.0)
        {
          result++;
        }
      }
      break;

    case 2:
      for (i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for (i[1] = 0; i[1] < A->N[1]; i[1]++)
        {
          if (spamm_get(i, A) != 0.0)
          {
            result++;
          }
        }
      }
      break;

    case 3:
      for (i[0] = 0; i[0] < A->N[0]; i[0]++) {
        for (i[1] = 0; i[1] < A->N[1]; i[1]++) {
          for (i[2] = 0; i[2] < A->N[2]; i[2]++)
          {
            if (spamm_get(i, A) != 0.0)
            {
              result++;
            }
          }
        }
      }
      break;

    default:
      SPAMM_FATAL("not implemented\n");
  }

  free(i);

  return result;
}

/** Print out some information on the matrix.
 *
 * @param A The matrix.
 */
void
spamm_print_info (const struct spamm_matrix_t *const A)
{
  int dim;
  unsigned int i;
  unsigned int N_matrix;
  unsigned int N_contiguous;

  assert(A != NULL);

  for (dim = 0, N_matrix = 1; dim < A->number_dimensions; dim++)
  {
    N_matrix *= A->N[dim];
  }

  for (i = 0, N_contiguous = A->N_padded; i < A->chunk_tier; i++)
  {
    N_contiguous >>= 1;
  }

  printf("number_dimensions = %u", A->number_dimensions);
  printf(", N = {");
  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    printf(" %u", A->N[dim]);
    if (dim+1 < A->number_dimensions)
    {
      printf(",");
    }
  }
  printf(" }");
  printf(", N_padded = %u", A->N_padded);
  printf(", chunk_tier = %u", A->chunk_tier);
  printf(", N_contiguous = %u", N_contiguous);
  printf(", use_linear_tree = %u", A->use_linear_tree);
  printf(", nnzero = %u (%1.2f %%)", spamm_number_nonzero(A), 100*(float) spamm_number_nonzero(A)/(float) N_matrix);
  printf("\n");
}
