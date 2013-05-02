/** @file
 *
 * Generate matrices for the tests.
 */

#include "test.h"

#include <spamm.h>
#include <stdlib.h>

/** Generate a random test matrix.
 *
 * @param number_dimensions The number of dimensions of the matrix.
 * @param is_sparse Whether the matrix is sparse or not.
 *
 * @return The newly allocated matrix.
 */
float *
generate_matrix (const unsigned int number_dimensions,
    const short is_sparse,
    unsigned int **const N)
{
  float *A;

  unsigned int dim;
  unsigned int *i;
  unsigned int N_contiguous;

  i = calloc(number_dimensions, sizeof(unsigned int));
  *N = calloc(number_dimensions, sizeof(unsigned int));

  for(dim = 0; dim < number_dimensions; dim++)
  {
    (*N)[dim] = 150+(int) ((0.5-(float) rand()/(float) RAND_MAX)*30);
  }

  N_contiguous = 1;
  for(dim = 0; dim < number_dimensions; dim++)
  {
    N_contiguous *= (*N)[dim];
  }

  A = (float*) calloc(N_contiguous, sizeof(float));

  if(is_sparse)
  {
    switch(number_dimensions)
    {
      case 1:
        for(i[0] = 0; i[0] < (*N)[0]; i[0]++)
        {
          A[i[0]] = rand()/(float) RAND_MAX;
        }
        break;

      case 2:
        for(i[0] = 0; i[0] < (*N)[0]; i[0]++) {
          for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < (*N)[1]; i[1]++)
          {
            A[spamm_index_row_major(i[0], i[1], (*N)[0], (*N)[1])] = rand()/(float) RAND_MAX;
          }
        }
        break;

      case 3:
        for(i[0] = 0; i[0] < (*N)[0]; i[0]++) {
          for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < (*N)[1]; i[1]++) {
            for((i[0] >= 10 ? i[2] = i[0]-10 : 0); i[2] < i[0]+10 && i[2] < (*N)[2]; i[2]++)
            {
              A[i[0]+(*N)[0]*(i[1]+(*N)[1]*i[2])] = rand()/(float) RAND_MAX;
            }
          }
        }
        break;

      default:
        SPAMM_FATAL("FIXME\n");
        break;
    }
  }

  else
  {
    for(i[0] = 0; i[0] < N_contiguous; i[0]++)
    {
      A[i[0]] = rand()/(float) RAND_MAX;
    }
  }

  free(i);

  return A;
}

/** Create a SpAMM matrix from a dense matrix.
 *
 * @param number_dimensions The number of dimensions.
 * @param N The size of each dimensions.
 *
 * @return The newly allocated SpAMM matrix.
 */
struct spamm_matrix_t *
create_spamm_from_dense (const unsigned int number_dimensions,
    const unsigned int *const N,
    const unsigned int chunk_tier,
    const short use_linear_tree,
    const float *const A_dense)
{
  struct spamm_matrix_t *A;
  unsigned int *i;

  A = spamm_new(number_dimensions, N, chunk_tier, use_linear_tree);
  i = calloc(number_dimensions, sizeof(unsigned int));

  switch(number_dimensions)
  {
    case 1:
      for(i[0] = 0; i[0] < N[0]; i[0]++)
      {
        spamm_set(i, A_dense[i[0]], A);
      }
      break;

    case 2:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++)
        {
          spamm_set(i, A_dense[i[0]*N[1]+i[1]], A);
        }
      }
      break;

    case 3:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++) {
          for(i[2] = 0; i[2] < N[2]; i[2]++)
          {
            spamm_set(i, A_dense[i[0]+N[0]*(i[1]+N[1]*i[2])], A);
          }
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

/** Compare a SpAMM matrix with a dense matrix.
 *
 * @return SPAMM_OK in case the matrices are identical, and SPAMM_ERROR in
 * case they are not.
 */
int
compare_spamm_to_dense (const struct spamm_matrix_t *const A,
    const float *const A_dense)
{
  int result = SPAMM_OK;

  unsigned int *i;
  unsigned int *N;

  i = calloc(spamm_get_number_dimensions(A), sizeof(unsigned int));
  N = spamm_get_N(A);

  switch(spamm_get_number_dimensions(A))
  {
    case 1:
      for(i[0] = 0; i[0] < N[0]; i[0]++)
      {
        if(A_dense[i[0]] != spamm_get(i, A))
        {
          result = SPAMM_ERROR;
          SPAMM_WARN("failed test at A[%u] (found %f, should be %f)\n",
              i[0], spamm_get(i, A), A_dense[i[0]]);
          break;
        }
        if(result) { break; }
      }
      break;

    case 2:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++)
        {
          if(A_dense[i[0]*N[1]+i[1]] != spamm_get(i, A))
          {
            result = SPAMM_ERROR;
            SPAMM_WARN("failed test at A[%u][%u] (found %f, should be %f)\n", i[0], i[1],
                spamm_get(i, A), A_dense[i[0]*N[1]+i[1]]);
            break;
          }
        }
        if(result) { break; }
      }
      break;

    case 3:
      for(i[0] = 0; i[0] < N[0]; i[0]++) {
        for(i[1] = 0; i[1] < N[1]; i[1]++) {
          for(i[2] = 0; i[2] < N[2]; i[2]++)
          {
            if(A_dense[i[0]+N[0]*(i[1]+N[1]*i[2])] != spamm_get(i, A))
            {
              result = SPAMM_ERROR;
              SPAMM_WARN("failed test at A[%u][%u][%u] (found %f, should be %f)\n", i[0], i[1], i[2],
                  spamm_get(i, A), A_dense[i[0]+N[0]*(i[1]+N[1]*i[2])]);
              break;
            }
          }
          if(result) { break; }
        }
        if(result) { break; }
      }
      break;

    default:
      SPAMM_FATAL("FIXME\n");
      break;
  }

  free(i);

  return SPAMM_OK;
}
