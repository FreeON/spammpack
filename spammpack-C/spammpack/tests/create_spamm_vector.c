#include <spamm.h>
#include <stdio.h>
#include <stdlib.h>

//#define PRINT_DEBUG

int
main (int argc, char **argv)
{
  int result = 0;
  unsigned int i[] = { 0 };
  const unsigned int N[] = { 511 };
  const unsigned int contiguous_tier = 5;
  const unsigned int N_block = 4;
  const short use_linear_tree = 1;
  struct spamm_matrix_t *A;
  float *A_dense = (float*) calloc(N[0], sizeof(float));

  enum spamm_layout_t layout = row_major;

  if (argc == 2)
  {
    layout = spamm_kernel_get_layout(argv[1]);
  }

  for (i[0] = 0; i[0] < N[0]; i[0]++)
  {
    A_dense[i[0]] = rand()/(float) RAND_MAX;
  }

  A = spamm_new(1, N, contiguous_tier, N_block, use_linear_tree);
  printf("A info: ");
  spamm_print_info(A);

  for (i[0] = 0; i[0] < N[0]; i[0]++)
  {
    spamm_set(i, A_dense[i[0]], A);
  }

#ifdef PRINT_DEBUG
  printf("A_dense:\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++)
  {
    printf(" %1.2f", A_dense[i[0]]);
  }
  printf("\n");

  /* For debugging, print out the whole tree. */
  printf("A:\n");
  spamm_print(A);
#endif

  for (i[0] = 0; i[0] < N[0]; i[0]++)
  {
    if (A_dense[i[0]] != spamm_get(i, A))
    {
      result = -1;
      printf("failed test at A[%u] (found %f, should be %f)\n", i[0],
          spamm_get(i, A), A_dense[i[0]]);
      break;
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
