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

  unsigned int N[] = { 13, 13 };

  const unsigned int linear_tier = 0;
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

  unsigned int max_i[] = { 0, 0 };
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

  spamm_print_info(A);
  spamm_print_info(B);
  spamm_print_info(C);

  //spamm_check(A, 1e-7);
  //spamm_check(B, 1e-7);
  //spamm_check(C, 1e-7);

#ifdef PRINT_DEBUG
  printf("A_dense =\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      printf(" %5.1f", A_dense[i[0]*N[1]+i[1]]);
    }
    printf("\n");
  }

  printf("A =\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      printf(" %5.1f", spamm_get(i, A));
    }
    printf("\n");
  }

  printf("B_dense =\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      printf(" %5.1f", B_dense[i[0]*N[1]+i[1]]);
    }
    printf("\n");
  }

  printf("C_dense =\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      printf(" %5.1f", C_dense[i[0]*N[1]+i[1]]);
    }
    printf("\n");
  }
#endif

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      C_dense[i[0]*N[1]+i[1]] *= beta;
    }
  }

  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++) {
      for (k = 0; k < N[0]; k++)
      {
        C_dense[i[0]*N[1]+i[1]] += alpha*A_dense[i[0]*N[1]+k]*B_dense[k*N[1]+i[1]];
      }
    }
  }

#ifdef PRINT_DEBUG
  printf("C_dense =\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      printf(" %5.1f", C_dense[i[0]*N[1]+i[1]]);
    }
    printf("\n");
  }
#endif

  timer = spamm_timer_new();
  spamm_timer_add_event(0x8000003b, timer);
  spamm_multiply(tolerance, alpha, A, B, beta, C, timer, NULL, kernel, NULL);
  spamm_timer_delete(&timer);

#ifdef PRINT_DEBUG
  printf("C =\n");
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++)
    {
      printf(" %5.1f", spamm_get(i, C));
    }
    printf("\n");
  }
#endif

#ifdef VERIFY_RESULT
  max_diff = 0;
  max_rel_diff = 0;
  for (i[0] = 0; i[0] < N[0]; i[0]++) {
    for (i[1] = 0; i[1] < N[1]; i[1]++) {
      for (k = 0; k < N[0]; k++)
      {
        if (fabs(C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C)) > max_diff)
        {
          max_diff = fabs(C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C));
          max_i[0] = i[0];
          max_i[1] = i[1];
        }

        if (C_dense[i[0]*N[1]+i[1]] != 0)
        {
          if (fabs((C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C))/C_dense[i[0]*N[1]+i[1]]) > max_rel_diff)
          {
            max_rel_diff = fabs((C_dense[i[0]*N[1]+i[1]]-spamm_get(i, C))/C_dense[i[0]*N[1]+i[1]]);
          }
        }
      }
    }
  }

  printf("max diff = %e, rel. diff = %e, A[%u][%u] = %e, A_reference[%u][%u] = %e\n",
      max_diff,
      (C_dense[max_i[0]*N[1]+max_i[1]] != 0.0 ? max_diff/C_dense[max_i[0]*N[1]+max_i[1]] : 0.0),
      max_i[0], max_i[1], spamm_get(max_i, C),
      max_i[0], max_i[1], C_dense[max_i[0]*N[1]+max_i[1]]);

  if (max_rel_diff > TEST_TOLERANCE)
  {
    printf("test failed\n");
    result = -1;
  }
#endif

  free(A_dense);
  free(B_dense);
  free(C_dense);

  spamm_delete(&A);
  spamm_delete(&B);
  spamm_delete(&C);

  return result;
}
