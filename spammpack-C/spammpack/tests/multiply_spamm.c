#include <spamm.h>
#include <math.h>
#include <stdio.h>

//#define RANDOM_ELEMENTS
#define PRINT_DEBUG

int
main ()
{
  int result = 0;

  unsigned int i, j, k;

  unsigned int N = 4;

  float alpha = 1.0;
  float beta = 1.0;

  float tolerance = 0.0;

  float *A_dense;
  float *B_dense;
  float *C_dense;

  struct spamm_t *A;
  struct spamm_t *B;
  struct spamm_t *C;

  unsigned int max_i, max_j;
  float max_diff;

  A_dense = (float*) malloc(sizeof(float)*N*N);
  B_dense = (float*) malloc(sizeof(float)*N*N);
  C_dense = (float*) malloc(sizeof(float)*N*N);

  A = spamm_new(N, N);
  B = spamm_new(N, N);
  C = spamm_new(N, N);

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
#ifdef RANDOM_ELEMENTS
      A_dense[i*N+j] = rand()/(float) RAND_MAX;
      B_dense[i*N+j] = rand()/(float) RAND_MAX;
      C_dense[i*N+j] = rand()/(float) RAND_MAX;
#else
      A_dense[i*N+j] = i*N+j;
      B_dense[i*N+j] = i*N+j;
      C_dense[i*N+j] = i*N+j;
#endif

      spamm_set(i, j, A_dense[i*N+j], A);
      spamm_set(i, j, B_dense[i*N+j], B);
      spamm_set(i, j, C_dense[i*N+j], C);
    }
  }

#ifdef PRINT_DEBUG
  printf("A_dense =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %5.1f", A_dense[i*N+j]);
    }
    printf("\n");
  }

  printf("B_dense =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %5.1f", B_dense[i*N+j]);
    }
    printf("\n");
  }

  printf("C_dense =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %5.1f", C_dense[i*N+j]);
    }
    printf("\n");
  }
#endif

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      C_dense[i*N+j] *= beta;
    }
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++)
      {
        C_dense[i*N+j] += alpha*A_dense[i*N+k]*B_dense[k*N+j];
      }
    }
  }

#ifdef PRINT_DEBUG
  printf("C_dense =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %5.1f", C_dense[i*N+j]);
    }
    printf("\n");
  }
#endif

  spamm_multiply(tolerance, alpha, A, B, beta, C);

#ifdef PRINT_DEBUG
  printf("C =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %5.1f", spamm_get(i, j, C));
    }
    printf("\n");
  }
#endif

  max_diff = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++)
      {
        if (fabs(C_dense[i*N+j]-spamm_get(i, j, C)) > max_diff)
        {
          max_diff = fabs(C_dense[i*N+j]-spamm_get(i, j, C));
          max_i = i;
          max_j = j;
        }
      }
    }
  }

  if (max_diff > 0)
  {
    printf("failed, max diff = %e for A(%u,%u)\n", max_diff, max_i, max_j);
  }

  free(A_dense);
  free(B_dense);
  free(C_dense);

  spamm_delete(&A);
  spamm_delete(&B);
  spamm_delete(&C);

  return result;
}
