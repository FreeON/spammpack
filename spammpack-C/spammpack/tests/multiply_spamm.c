#include <spamm.h>
#include <math.h>

int
main ()
{
  int result = 0;

  unsigned int i, j, k;

  unsigned int N = 64;

  float alpha = 1.2;
  float beta = 0.5;

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
      A_dense[i*N+j] = rand()/(float) RAND_MAX;
      B_dense[i*N+j] = rand()/(float) RAND_MAX;
      C_dense[i*N+j] = rand()/(float) RAND_MAX;

      spamm_set(i, j, A_dense[i*N+j], A);
      spamm_set(i, j, B_dense[i*N+j], B);
      spamm_set(i, j, C_dense[i*N+j], C);
    }
  }

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
        C_dense[i*N+j] += alpha*A_dense[i*N+j]*B_dense[i*N+j];
      }
    }
  }

  spamm_multiply(alpha, A, B, beta, C);

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
