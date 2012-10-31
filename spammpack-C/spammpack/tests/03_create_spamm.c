#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

int
main (int argc, char **argv)
{
  int result = SPAMM_OK;

  unsigned int number_dimensions;

  unsigned int *i;
  unsigned int *N;

  unsigned int N_contiguous;

  unsigned int dim;

  const unsigned int contiguous_tier = 4;
  const unsigned int N_block = 4;

  short use_linear_tree;
  short is_sparse;
  struct spamm_matrix_t *A;

  float *A_dense;

  for (number_dimensions = 1; number_dimensions < 3; number_dimensions++)
  {
    i = calloc(number_dimensions, sizeof(unsigned int));
    N = calloc(number_dimensions, sizeof(unsigned int));

    for (dim = 0; dim < number_dimensions; dim++)
    {
      N[dim] = 512+(int) ((0.5-(float) rand()/(float) RAND_MAX)*20);
    }

    N_contiguous = 1;
    for (dim = 0; dim < number_dimensions; dim++)
    {
      N_contiguous *= N[dim];
    }

    A_dense = (float*) calloc(N_contiguous, sizeof(float));

    for (use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
    {
      for (is_sparse = 0; is_sparse < 2; is_sparse++)
      {
        printf("dim: %u, linTree: %u, sparse: %u, ", number_dimensions, use_linear_tree, is_sparse);
        if (is_sparse)
        {
          switch (number_dimensions)
          {
            case 1:
              for (i[0] = 0; i[0] < N[0]; i[0]++)
              {
                A_dense[i[0]] = rand()/(float) RAND_MAX;
              }
              break;

            case 2:
              for (i[0] = 0; i[0] < N[0]; i[0]++) {
                for ((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++)
                {
                  A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = rand()/(float) RAND_MAX;
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
          for (i[0] = 0; i[0] < N_contiguous; i[0]++)
          {
            A_dense[i[0]] = rand()/(float) RAND_MAX;
          }
        }

        A = spamm_new(number_dimensions, N, contiguous_tier, N_block, use_linear_tree);
        printf("A info: ");
        spamm_print_info(A);

        switch (number_dimensions)
        {
          case 1:
            for (i[0] = 0; i[0] < N[0]; i[0]++)
            {
              spamm_set(i, A_dense[i[0]], A);
            }
            break;

          case 2:
            for (i[0] = 0; i[0] < N[0]; i[0]++) {
              for (i[1] = 0; i[1] < N[1]; i[1]++)
              {
                spamm_set(i, A_dense[i[0]*N[1]+i[1]], A);
              }
            }

#ifdef PRINT_DEBUG
            printf("A_dense:\n");
            for (i[0] = 0; i[0] < N[0]; i[0]++) {
              for (i[1] = 0; i[1] < N[1]; i[1]++)
              {
                printf(" %1.2f", A_dense[i[0]*N[1]+i[1]]);
              }
              printf("\n");
            }

            /* For debugging, print out the whole tree. */
            printf("A:\n");
            spamm_print(A);
#endif
            break;

          default:
            SPAMM_FATAL("FIXME\n");
            break;
        }

        switch (number_dimensions)
        {
          case 1:
            for (i[0] = 0; i[0] < N[0]; i[0]++)
            {
              if (A_dense[i[0]] != spamm_get(i, A))
              {
                result = SPAMM_ERROR;
                printf("failed test at A[%u] (found %f, should be %f)\n",
                    i[0], spamm_get(i, A), A_dense[i[0]]);
                break;
              }
              if (result) { break; }
            }
            break;

          case 2:
            for (i[0] = 0; i[0] < N[0]; i[0]++) {
              for (i[1] = 0; i[1] < N[1]; i[1]++)
              {
                if (A_dense[i[0]*N[1]+i[1]] != spamm_get(i, A))
                {
                  result = SPAMM_ERROR;
                  printf("failed test at A[%u][%u] (found %f, should be %f)\n", i[0], i[1],
                      spamm_get(i, A), A_dense[i[0]*N[1]+i[1]]);
                  break;
                }
              }
              if (result) { break; }
            }
            break;

          default:
            SPAMM_FATAL("FIXME\n");
            break;
        }

#ifdef PRINT_DEBUG
        printf("test passed\n");
#endif

        spamm_delete(&A);
      }
    }
    free(i);
    free(N);
    free(A_dense);
  }

  return result;
}
