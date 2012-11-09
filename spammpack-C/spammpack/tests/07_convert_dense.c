#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

int
main (int argc, char **argv)
{
  int result = 0;
  unsigned int i[2];
  const unsigned int N[] = { 1024, 1024 };
  const unsigned int chunk_tier = 6;

  const short use_linear_tree = 1;
  struct spamm_matrix_t *A = NULL;
  float *A_dense = (float*) malloc(sizeof(float)*N[0]*N[1]);

  enum spamm_layout_t layout = row_major;

  if (argc == 2)
  {
    layout = spamm_kernel_get_layout(argv[1]);
  }

  for (i[0] = 0; i[0] < N[0]*N[1]; i[0]++)
  {
    A_dense[i[0]] = rand()/(float) RAND_MAX;
  }

  A = spamm_convert_dense_to_spamm(2, N, chunk_tier, use_linear_tree, row_major, A_dense);

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
  spamm_print(A);
#endif

  /* Check matrix. */
  if ((result = spamm_check(A, 1e-7)) != SPAMM_OK)
  {
    printf("failed spamm_hashed_check()\n");
    return result;
  }

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      if (A_dense[i[0]*N[1]+i[1]] != spamm_get(i, A))
      {
        result = -1;
        printf("failed test at A[%u][%u] (found %f, should be %f)\n", i[0], i[1],
            spamm_get(i, A), A_dense[i[0]*N[1]+i[0]]);
        break;
      }
    }
    if (result < 0) { break; }
  }

#ifdef PRINT_DEBUG
  printf("test passed\n");
#endif

  spamm_delete(&A);
  free(A_dense);
  return result;
}
