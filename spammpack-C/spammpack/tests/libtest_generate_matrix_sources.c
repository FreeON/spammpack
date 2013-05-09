/** @file */

#include "test.h"
#include <spamm.h>

#include <math.h>
#include <stdlib.h>

#define SPARSITY 0.9

/** Generate a random test matrix. The matrix elements are stored in row-major
 * order.
 *
 * @param number_dimensions The number of dimensions of the matrix.
 * @param matrix_type The matrix type.
 * @param N The shape of the matrix.
 * @param gamma The decay factor for exponential-decay type matrices.
 *
 * @return The newly allocated matrix.
 */
MATRIX_TYPE *
SPAMM_FUNCTION(generate_matrix, MATRIX_TYPE) (const unsigned int number_dimensions,
    const enum matrix_t matrix_type,
    const unsigned int *const N,
    const MATRIX_TYPE gamma)
{
  MATRIX_TYPE *A;

  int *i;
  int dim;

  short is_square;

  unsigned int N_linear;
  unsigned int nonzero_elements;

  for(dim = 1, is_square = 1; dim < number_dimensions; dim++)
  {
    if(N[dim-1] != N[dim])
    {
      is_square = 0;
      break;
    }
  }

  N_linear = 1;
  for(dim = 0; dim < number_dimensions; dim++)
  {
    N_linear *= N[dim];
  }

  i = calloc(number_dimensions, sizeof(int));
  A = (MATRIX_TYPE*) calloc(N_linear, sizeof(MATRIX_TYPE));

  switch(matrix_type)
  {
    case full:
      for(i[0] = 0; i[0] < N_linear; i[0]++)
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
              A[i[1]+N[1]*i[0]] = rand()/(MATRIX_TYPE) RAND_MAX;
            }
          }
          break;

        case 3:
          for(i[0] = 0; i[0] < N[0]; i[0]++) {
            for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++) {
              for((i[0] >= 10 ? i[2] = i[0]-10 : 0); i[2] < i[0]+10 && i[2] < N[2]; i[2]++)
              {
                A[i[2]+N[2]*(i[1]+N[1]*i[0])] = rand()/(MATRIX_TYPE) RAND_MAX;
              }
            }
          }
          break;

        default:
          SPAMM_FATAL("FIXME\n");
          break;
      }
      break;

    case exponential_decay:
      if(!is_square)
      {
        SPAMM_FATAL("FIXME\n");
      }

      switch(number_dimensions)
      {
        case 1:
          for(i[0] = 0; i[0] < N[0]; i[0]++)
          {
            A[i[0]] = 0.8+rand()/(double) RAND_MAX;
          }
          break;

        case 2:
          for(i[0] = 0; i[0] < N[0]; i[0]++)
          {
            A[i[0]+N[1]*i[0]] = 0.8+rand()/(double) RAND_MAX;
          }

          for(i[0] = 0; i[0] < N[0]; i[0]++) {
            for(i[1] = i[0]+1; i[1] < N[1]; i[1]++)
            {
              /* Exponential decay. */
#if MATRIX_TYPE == float
              A[i[1]+N[1]*i[0]] = A[i[0]+N[1]*i[0]]*expf(-fabsf(i[0]-i[1])/gamma);
#else
              A[i[1]+N[1]*i[0]] = A[i[0]+N[1]*i[0]]*exp(-fabs(i[0]-i[1])/gamma);
#endif
              A[i[0]+N[1]*i[1]] = A[i[1]+N[1]*i[0]];
            }
          }
          break;

        case 3:
          for(i[0] = 0; i[0] < N[0]; i[0]++)
          {
            A[i[0]+N[2]*(i[0]+N[1]*i[0])] = 0.8+rand()/(double) RAND_MAX;
          }

          for(i[0] = 0; i[0] < N[0]; i[0]++) {
            for(i[1] = i[0]+1; i[1] < N[1]; i[1]++) {
              for(i[2] = i[1]+1; i[2] < N[2]; i[2]++)
              {
                /* Exponential decay. */
#if MATRIX_TYPE == float
                A[i[2]+N[2]*(i[1]+N[1]*i[0])] = A[i[0]+N[2]*(i[0]+N[1]*i[0])]*expf(-fabsf(i[0]-i[1])/gamma);
#else
                A[i[2]+N[2]*(i[1]+N[1]*i[0])] = A[i[0]+N[2]*(i[0]+N[1]*i[0])]*exp(-fabs(i[0]-i[1])/gamma);
#endif
                A[i[2]+N[2]*(i[0]+N[1]*i[1])] = A[i[2]+N[2]*(i[1]+N[1]*i[0])];
                A[i[1]+N[2]*(i[0]+N[1]*i[2])] = A[i[2]+N[2]*(i[1]+N[1]*i[0])];
                A[i[1]+N[2]*(i[2]+N[1]*i[0])] = A[i[2]+N[2]*(i[1]+N[1]*i[0])];
                A[i[0]+N[2]*(i[1]+N[1]*i[2])] = A[i[2]+N[2]*(i[1]+N[1]*i[0])];
                A[i[0]+N[2]*(i[2]+N[1]*i[1])] = A[i[2]+N[2]*(i[1]+N[1]*i[0])];
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
      for(nonzero_elements = 0; nonzero_elements < (unsigned int) ceil((1-SPARSITY)*N_linear); )
      {
        i[0] = (unsigned int) floor(rand()/(double) RAND_MAX * N_linear);
        if(A[i[0]] == 0.0)
        {
          A[i[0]] = rand()/(MATRIX_TYPE) RAND_MAX;
          nonzero_elements++;
        }
      }
      break;

    default:
      SPAMM_FATAL("FIXME\n");
      break;
  }
  free(i);

  return A;
}
