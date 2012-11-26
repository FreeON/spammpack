#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

int
main (int argc, char **argv)
{
  int result = 0;

  float beta = 1.2;

  unsigned int number_dimensions;

  unsigned int *i;
  unsigned int *N;

  unsigned int N_contiguous;

  unsigned int dim;

  const unsigned int chunk_tier = 5;

  short use_linear_tree;
  short is_sparse;

  struct spamm_matrix_t *A = NULL;
  struct spamm_matrix_t *B = NULL;

  float *B_dense = (float*) calloc(N[0], sizeof(float));

  for(number_dimensions = 1; number_dimensions <= 3; number_dimensions++)
  {
    i = calloc(number_dimensions, sizeof(unsigned int));
    N = calloc(number_dimensions, sizeof(unsigned int));

    for(dim = 0; dim < number_dimensions; dim++)
    {
      N[dim] = 150+(int) ((0.5-(float) rand()/(float) RAND_MAX)*30);
    }

    N_contiguous = 1;
    for(dim = 0; dim < number_dimensions; dim++)
    {
      N_contiguous *= N[dim];
    }

    for(use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
    {
      for(is_sparse = 0; is_sparse < 2; is_sparse++)
      {
        B_dense = (float*) calloc(N_contiguous, sizeof(float));

        printf("dim: %u, linTree: %u, sparse: %u\n", number_dimensions, use_linear_tree, is_sparse);
        if(is_sparse)
        {
          switch(number_dimensions)
          {
            case 1:
              for(i[0] = 0; i[0] < N[0]; i[0]++)
              {
                B_dense[i[0]] = rand()/(float) RAND_MAX;
              }
              break;

            case 2:
              for(i[0] = 0; i[0] < N[0]; i[0]++) {
                for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++)
                {
                  B_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = rand()/(float) RAND_MAX;
                }
              }
              break;

            case 3:
              for(i[0] = 0; i[0] < N[0]; i[0]++) {
                for((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++) {
                  for((i[0] >= 10 ? i[2] = i[0]-10 : 0); i[2] < i[0]+10 && i[2] < N[2]; i[2]++)
                  {
                    B_dense[i[0]+N[0]*(i[1]+N[1]*i[2])] = rand()/(float) RAND_MAX;
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
            B_dense[i[0]] = rand()/(float) RAND_MAX;
          }
        }

        B = spamm_new(number_dimensions, N, chunk_tier, use_linear_tree);

        switch(number_dimensions)
        {
          case 1:
            for(i[0] = 0; i[0] < N[0]; i[0]++)
            {
              spamm_set(i, B_dense[i[0]], B);
            }
            break;

          case 2:
            for(i[0] = 0; i[0] < N[0]; i[0]++) {
              for(i[1] = 0; i[1] < N[1]; i[1]++)
              {
                spamm_set(i, B_dense[i[0]*N[1]+i[1]], B);
              }
            }

#ifdef PRINT_DEBUG
            printf("B_dense:\n");
            for(i[0] = 0; i[0] < N[0]; i[0]++) {
              for(i[1] = 0; i[1] < N[1]; i[1]++)
              {
                printf(" %1.2f", B_dense[i[0]*N[1]+i[1]]);
              }
              printf("\n");
            }

            /* For debugging, print out the whole tree. */
            printf("A:\n");
            spamm_print(A);
#endif
            break;

          case 3:
            for(i[0] = 0; i[0] < N[0]; i[0]++) {
              for(i[1] = 0; i[1] < N[1]; i[1]++) {
                for(i[2] = 0; i[2] < N[2]; i[2]++)
                {
                  spamm_set(i, B_dense[i[0]+N[0]*(i[1]+N[1]*i[2])], B);
                }
              }
            }
            break;

          default:
            SPAMM_FATAL("FIXME\n");
            break;
        }

        printf("B info: ");
        spamm_print_info(B);
        spamm_check(B, 1e-8);

        spamm_copy(&A, beta, B);

        printf("A info: ");
        spamm_print_info(A);
        spamm_check(A, 1e-6);

        switch(number_dimensions)
        {
          case 1:
            for(i[0] = 0; i[0] < N[0]; i[0]++)
            {
              if(beta*B_dense[i[0]] != spamm_get(i, A))
              {
                result = -1;
                printf("failed test at A[%u] (found %f, should be %f)\n", i[0],
                    spamm_get(i, A), beta*B_dense[i[0]]);
                break;
              }
              if(result < 0) { break; }
            }
            break;

          case 2:
            for(i[0] = 0; i[0] < N[0]; i[0]++) {
              for(i[1] = 0; i[1] < N[1]; i[1]++)
              {
                if(beta*B_dense[i[0]*N[1]+i[1]] != spamm_get(i, A))
                {
                  result = -1;
                  printf("failed test at A[%u][%u] (found %f, should be %f)\n", i[0], i[1],
                      spamm_get(i, A), beta*B_dense[i[0]*N[1]+i[1]]);
                  break;
                }
                if(result < 0) { break; }
              }
            }
            break;

          case 3:
            for(i[0] = 0; i[0] < N[0]; i[0]++) {
              for(i[1] = 0; i[1] < N[1]; i[1]++) {
                for(i[2] = 0; i[2] < N[2]; i[2]++)
                {
                  if(beta*B_dense[i[0]+N[0]*(i[1]+N[1]*i[2])] != spamm_get(i, A))
                  {
                    result = -1;
                    printf("failed test at A[%u][%u][%u] (found %f, should be %f)\n", i[0], i[1], i[2],
                        spamm_get(i, A), beta*B_dense[i[0]+N[0]*(i[1]+N[1]*i[2])]);
                    break;
                  }
                  if(result < 0) { break; }
                }
              }
            }
            break;

          default:
            SPAMM_FATAL("FIXME\n");
            break;
        }

#ifdef PRINT_DEBUG
        printf("test passed\n");
#endif

        free(B_dense);
        spamm_delete(&A);
        spamm_delete(&B);
      }
    }
  }

  return result;
}
