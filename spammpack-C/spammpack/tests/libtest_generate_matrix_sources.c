/** @file */

#include "test.h"
#include "spamm.h"

#include <stdlib.h>

/** Generate a random test matrix.
 *
 * @param number_dimensions The number of dimensions of the matrix.
 * @param matrix_type The matrix type.
 * @param N The shape of the matrix.
 *
 * @return The newly allocated matrix.
 */
MATRIX_TYPE *
SPAMM_FUNCTION(generate_matrix, MATRIX_TYPE) (const unsigned int number_dimensions,
    const enum matrix_t matrix_type,
    const unsigned int *const N)
{
  MATRIX_TYPE *A;

  unsigned int dim;
  unsigned int *i;
  unsigned int N_contiguous;

  N_contiguous = 1;
  for(dim = 0; dim < number_dimensions; dim++)
  {
    N_contiguous *= N[dim];
  }

  i = calloc(number_dimensions, sizeof(unsigned int));
  A = (MATRIX_TYPE*) calloc(N_contiguous, sizeof(MATRIX_TYPE));

  switch(matrix_type)
  {
    case full:
      for(i[0] = 0; i[0] < N_contiguous; i[0]++)
      {
        A[i[0]] = rand()/(MATRIX_TYPE) RAND_MAX;
      }
      break;

    case diagonally_banded:
      switch(number_dimensions)
      {
        case 1:
          for(i[0] = 0; i[0] < N[0]; i[0]++)
          {
            A[i[0]] = rand()/(MATRIX_TYPE) RAND_MAX;
          }
          break;

        case 2:
          for(i[0] = 0; i[0] < N[0]; i[0]++) {
            for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++)
            {
              A[spamm_index_row_major(i[0], i[1], N[0], N[1])] = rand()/(MATRIX_TYPE) RAND_MAX;
            }
          }
          break;

        case 3:
          for(i[0] = 0; i[0] < N[0]; i[0]++) {
            for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++) {
              for((i[0] >= 10 ? i[2] = i[0]-10 : 0); i[2] < i[0]+10 && i[2] < N[2]; i[2]++)
              {
                A[i[0]+N[0]*(i[1]+N[1]*i[2])] = rand()/(MATRIX_TYPE) RAND_MAX;
              }
            }
          }
          break;

        default:
          SPAMM_FATAL("FIXME\n");
          break;
      }
      break;

    case sparse_random:
      break;

    default:
      SPAMM_FATAL("FIXME\n");
      break;
  }
  free(i);

  return A;
}
