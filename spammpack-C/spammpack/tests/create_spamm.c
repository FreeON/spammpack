#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

int
main (int argc, char **argv)
{
  int result = 0;
  unsigned int i, j;
  const unsigned int N = 78;
  const unsigned int linear_tier = 1;
  const unsigned int contiguous_tier = 4;
  struct spamm_matrix_t *A;
  float *A_dense = (float*) calloc(N*N, sizeof(float));

  enum spamm_layout_t layout = row_major;

  if (argc == 2)
  {
    layout = spamm_kernel_get_layout(argv[1]);
  }

  for (i = 0; i < N*N; i++)
  {
    A_dense[i] = rand()/(float) RAND_MAX;
  }

  A = spamm_new(N, N, linear_tier, contiguous_tier, layout);
  printf("A info: ");
  spamm_print_info(A);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      spamm_set(i, j, A_dense[i*N+j], A);
    }
  }

#ifdef PRINT_DEBUG
  printf("A_dense:\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %1.2f", A_dense[i*N+j]);
    }
    printf("\n");
  }

  /* For debugging, print out the whole tree. */
  printf("A:\n");
  spamm_print(A);
#endif

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (A_dense[i*N+j] != spamm_get(i, j, A))
      {
        result = -1;
        printf("failed test at A[%u][%u] (found %f, should be %f)\n", i, j,
            spamm_get(i, j, A), A_dense[i*N+j]);
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
