#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

int
main (int argc, char **argv)
{
  int result = 0;
  unsigned int i[2];
  const unsigned int N[] = { 511, 513 };
  const unsigned int contiguous_tier = 4;
  const unsigned int N_block = 4;
  short use_linear_tree;
  short is_sparse;
  struct spamm_matrix_t *A;

  float *A_dense = (float*) calloc(N[0]*N[1], sizeof(float));

  for (use_linear_tree = 0; use_linear_tree < 2; use_linear_tree++)
  {
    for (is_sparse = 0; is_sparse < 2; is_sparse++)
    {
      printf("use_linear_tree: %u, is_sparse: %u, ", use_linear_tree, is_sparse);
      if (is_sparse)
      {
        for (i[0] = 0; i[0] < N[0]; i[0]++)
        {
          for ((i[0] >= 10 ? i[1] = i[0]-10 : 0); i[1] < i[0]+10 && i[1] < N[1]; i[1]++)
          {
            A_dense[spamm_index_row_major(i[0], i[1], N[0], N[1])] = rand()/(float) RAND_MAX;
          }
        }
      }

      else
      {
        for (i[0] = 0; i[0] < N[0]*N[1]; i[0]++)
        {
          A_dense[i[0]] = rand()/(float) RAND_MAX;
        }
      }

      A = spamm_new(2, N, contiguous_tier, N_block, use_linear_tree);
      printf("A info: ");
      spamm_print_info(A);

      for (i[0] = 0; i[0] < N[0]; i[0]++) {
        for (i[1] = 0; i[1] < N[1]; i[1]++)
        {
          //printf("%u %u %f\n", i[0], i[1], A_dense[i[0]*N[1]+i[1]]);
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

      for (i[0] = 0; i[0] < N[0]; i[0]++) {
        for (i[1] = 0; i[1] < N[1]; i[1]++)
        {
          if (A_dense[i[0]*N[1]+i[1]] != spamm_get(i, A))
          {
            result = -1;
            printf("failed test at A[%u][%u] (found %f, should be %f)\n", i[0], i[1],
                spamm_get(i, A), A_dense[i[0]*N[1]+i[1]]);
            break;
          }
        }
        if (result < 0) { break; }
      }

#ifdef PRINT_DEBUG
      printf("test passed\n");
#endif

      spamm_delete(&A);
      break;
    }
    break;
  }
  free(A_dense);

  return result;
}
