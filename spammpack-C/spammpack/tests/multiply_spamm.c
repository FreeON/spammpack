#include <spamm.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define VERIFY_RESULT
#define RANDOM_ELEMENTS
//#define PRINT_DEBUG

#define TEST_TOLERANCE 5e-7

int
main (int argc, char **argv)
{
  int result = 0;

  unsigned int i[2];
  unsigned int k;

  unsigned int N[] = { 513, 513 };

  const unsigned int linear_tier = 4;
  const unsigned int contiguous_tier = 5;

  double alpha = 1.2;
  double beta = 0.5;

  float tolerance = 0.0;

  double *A_dense;
  double *B_dense;
  double *C_dense;

  struct spamm_matrix_t *A;
  struct spamm_matrix_t *B;
  struct spamm_matrix_t *C;

  unsigned int max_i = 0;
  unsigned int max_j = 0;
  double max_diff;
  double max_rel_diff;

  enum spamm_kernel_t kernel = kernel_standard_SSE;
  struct spamm_timer_t *timer;

  if (argc == 2)
  {
    kernel = spamm_kernel_get_kernel(argv[1]);
  }

  A_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  B_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);
  C_dense = (double*) malloc(sizeof(double)*N[0]*N[1]);

  A = spamm_new(2, N, linear_tier, contiguous_tier, spamm_kernel_suggest_layout(kernel));
  B = spamm_new(2, N, linear_tier, contiguous_tier, spamm_kernel_suggest_layout(kernel));
  C = spamm_new(2, N, linear_tier, contiguous_tier, spamm_kernel_suggest_layout(kernel));

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
#ifdef RANDOM_ELEMENTS
      A_dense[i[0]*N[1]+i[1]] = rand()/(double) RAND_MAX;
      B_dense[i[0]*N[1]+i[1]] = rand()/(double) RAND_MAX;
      C_dense[i[0]*N[1]+i[1]] = rand()/(double) RAND_MAX;
#else
      A_dense[i[0]*N[1]+i[1]] = i[0]*N[1]+i[1];
      B_dense[i[0]*N[1]+i[1]] = i[0]*N[1]+i[1];
      C_dense[i[0]*N[1]+i[1]] = i[0]*N[1]+i[1];
#endif

      spamm_set(i, A_dense[i[0]*N[1]+i[1]], A);
      spamm_set(i, B_dense[i[0]*N[1]+i[1]], B);
      spamm_set(i, C_dense[i[0]*N[1]+i[1]], C);
    }
  }

//j  //spamm_check(A, 1e-7);
//j  //spamm_check(B, 1e-7);
//j  //spamm_check(C, 1e-7);
//j
//j#ifdef PRINT_DEBUG
//j  printf("A_dense =\n");
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      printf(" %5.1f", A_dense[i*N+j]);
//j    }
//j    printf("\n");
//j  }
//j
//j  printf("A =\n");
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      printf(" %5.1f", spamm_get(i, j, A));
//j    }
//j    printf("\n");
//j  }
//j
//j  printf("B_dense =\n");
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      printf(" %5.1f", B_dense[i*N+j]);
//j    }
//j    printf("\n");
//j  }
//j
//j  printf("C_dense =\n");
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      printf(" %5.1f", C_dense[i*N+j]);
//j    }
//j    printf("\n");
//j  }
//j#endif
//j
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      C_dense[i*N+j] *= beta;
//j    }
//j  }
//j
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++) {
//j      for (k = 0; k < N; k++)
//j      {
//j        C_dense[i*N+j] += alpha*A_dense[i*N+k]*B_dense[k*N+j];
//j      }
//j    }
//j  }
//j
//j#ifdef PRINT_DEBUG
//j  printf("C_dense =\n");
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      printf(" %5.1f", C_dense[i*N+j]);
//j    }
//j    printf("\n");
//j  }
//j#endif
//j
//j  timer = spamm_timer_new();
//j  spamm_timer_add_event(0x8000003b, timer);
//j  spamm_matrix_multiply(tolerance, alpha, A, B, beta, C, timer, kernel);
//j  spamm_timer_delete(&timer);
//j
//j#ifdef PRINT_DEBUG
//j  printf("C =\n");
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++)
//j    {
//j      printf(" %5.1f", spamm_get(i, j, C));
//j    }
//j    printf("\n");
//j  }
//j#endif
//j
//j#ifdef VERIFY_RESULT
//j  max_diff = 0;
//j  max_rel_diff = 0;
//j  for (i = 0; i < N; i++) {
//j    for (j = 0; j < N; j++) {
//j      for (k = 0; k < N; k++)
//j      {
//j        if (fabs(C_dense[i*N+j]-spamm_get(i, j, C)) > max_diff)
//j        {
//j          max_diff = fabs(C_dense[i*N+j]-spamm_get(i, j, C));
//j          max_i = i;
//j          max_j = j;
//j        }
//j
//j        if (C_dense[i*N+j] != 0)
//j        {
//j          if (fabs((C_dense[i*N+j]-spamm_get(i, j, C))/C_dense[i*N+j]) > max_rel_diff)
//j          {
//j            max_rel_diff = fabs((C_dense[i*N+j]-spamm_get(i, j, C))/C_dense[i*N+j]);
//j          }
//j        }
//j      }
//j    }
//j  }
//j
//j  printf("max diff = %e, rel. diff = %e, A(%u,%u) = %e, A_reference(%u,%u) = %e\n",
//j      max_diff,
//j      (C_dense[max_i*N+max_j] != 0.0 ? max_diff/C_dense[max_i*N+max_j] : 0.0),
//j      max_i, max_j, spamm_get(max_i, max_j, C),
//j      max_i, max_j, C_dense[max_i*N+max_j]);
//j
//j  if (max_rel_diff > TEST_TOLERANCE)
//j  {
//j    printf("test failed\n");
//j    result = -1;
//j  }
//j#endif
//j
//j  free(A_dense);
//j  free(B_dense);
//j  free(C_dense);
//j
//j  spamm_matrix_delete(&A);
//j  spamm_matrix_delete(&B);
//j  spamm_matrix_delete(&C);

  return result;
}
