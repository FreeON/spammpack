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

    default:
      SPAMM_FATAL("not implemented\n");
  }

  free(i);

  return result;
}

/** Return the total memory consumption of a matrix.
 *
 * @param A The matrix.
 *
 * @return The number of bytes allocated for this matrix.
 */
unsigned int
spamm_memory (const struct spamm_hashed_t *A)
{
  unsigned int total = 0;

  total = sizeof(struct spamm_hashed_t);
  total += sizeof(struct spamm_hashtable_t*)*(A->kernel_tier+1);

  return total;
}

/** Print out some information on the matrix.
 *
 * @param A The matrix.
 */
void
spamm_print_info (const struct spamm_matrix_t *const A)
{
  int dim;

  assert(A != NULL);

  printf("number_dimensions = %u", A->number_dimensions);
  for (dim = 0; dim < A->number_dimensions; dim++)
  {
    printf(", N[%i] = %u", dim, A->N[dim]);
  }
  printf(", N_padded = %u", A->N_padded);
  printf(", depth = %u", A->depth);
  printf(", contiguous_tier = %u", A->contiguous_tier);
  printf(", N_contiguous = %u", 0);
  printf(", kernel_tier = %u", A->kernel_tier);
  printf(", N_block = %u", A->N_block);
  printf(", use_linear_tree = %u", A->use_linear_tree);

  printf("\n");
}
