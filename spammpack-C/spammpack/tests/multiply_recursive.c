#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "config.h"
#include "spamm.h"

#define TEST_TOLERANCE 5e-7

#ifndef SGEMM
#define SGEMM sgemm_
#endif

//#define PRINT_DEBUG

int
main ()
{
  int result = 0;

  const int N =  513;
  const int blocksize = 16;
  const float alpha = 1.2;
  const float beta = 0.5;

  float tolerance = 1e-4;

  struct spamm_recursive_t *A = NULL;
  struct spamm_recursive_t *B = NULL;
  struct spamm_recursive_t *C = NULL;
  float *A_dense = NULL;
  float *B_dense = NULL;
  float *C_dense = NULL;

  int i, j, k;

  unsigned int max_i = 0;
  unsigned int max_j = 0;
  double max_diff;
  double max_rel_diff;

  unsigned int recursive_number_products = 0;

  A = spamm_recursive_new(N, N, blocksize);
  B = spamm_recursive_new(N, N, blocksize);
  C = spamm_recursive_new(N, N, blocksize);

  A_dense = calloc(N*N, sizeof(float));
  B_dense = calloc(N*N, sizeof(float));
  C_dense = calloc(N*N, sizeof(float));

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      A_dense[spamm_index_row_major(i, j, N, N)] = rand()/(float) RAND_MAX;
      B_dense[spamm_index_row_major(i, j, N, N)] = rand()/(float) RAND_MAX;
      C_dense[spamm_index_row_major(i, j, N, N)] = rand()/(float) RAND_MAX;

      spamm_recursive_set(i, j, A_dense[spamm_index_row_major(i, j, N, N)], A);
      spamm_recursive_set(i, j, B_dense[spamm_index_row_major(i, j, N, N)], B);
      spamm_recursive_set(i, j, C_dense[spamm_index_row_major(i, j, N, N)], C);
    }
  }


  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (A_dense[spamm_index_row_major(i, j, N, N)] != spamm_recursive_get(i, j, A))
      {
        printf("[%s:%i] element mismatch\n", __FILE__, __LINE__);
        return -1;
      }

      if (B_dense[spamm_index_row_major(i, j, N, N)] != spamm_recursive_get(i, j, B))
      {
        printf("[%s:%i] element mismatch\n", __FILE__, __LINE__);
        return -1;
      }

      if (C_dense[spamm_index_row_major(i, j, N, N)] != spamm_recursive_get(i, j, C))
      {
        printf("[%s:%i] element mismatch\n", __FILE__, __LINE__);
        return -1;
      }
    }
  }

#ifdef PRINT_DEBUG
  /* Print matrices for debugging. */
  printf("alpha = %1.2f\n", alpha);
  printf("beta = %1.2f\n", beta);
  printf("A =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", A_dense[spamm_index_row_major(i, j, N, N)]);
    }
    printf("\n");
  }

  printf("B =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", B_dense[spamm_index_row_major(i, j, N, N)]);
    }
    printf("\n");
  }

  printf("C(before) =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", C_dense[spamm_index_row_major(i, j, N, N)]);
    }
    printf("\n");
  }
#endif

  /* Multiply by beta. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      C_dense[spamm_index_row_major(i, j, N, N)] *= beta;
    }
  }

  /* Multiply A and B. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++)
      {
        C_dense[spamm_index_row_major(i, j, N, N)] += alpha*
          A_dense[spamm_index_row_major(i, k, N, N)]*
          B_dense[spamm_index_row_major(k, j, N, N)];
      }
    }
  }

#ifdef PRINT_DEBUG
  printf("C(after) =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", C_dense[spamm_index_row_major(i, j, N, N)]);
    }
    printf("\n");
  }
#endif

  spamm_recursive_multiply(tolerance, alpha, A, B, beta, C, NULL, SGEMM, &recursive_number_products);

#ifdef PRINT_DEBUG
  printf("C(SpAMM) =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", spamm_recursive_get(i, j, C));
    }
    printf("\n");
  }
#endif

  /* Compare results. */
  max_diff = 0;
  max_rel_diff = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++)
      {
        if (fabs(C_dense[i*N+j]-spamm_recursive_get(i, j, C)) > max_diff)
        {
          max_diff = fabs(C_dense[i*N+j]-spamm_recursive_get(i, j, C));
          max_i = i;
          max_j = j;
        }

        if (C_dense[i*N+j] != 0)
        {
          if (fabs((C_dense[i*N+j]-spamm_recursive_get(i, j, C))/C_dense[i*N+j]) > max_rel_diff)
          {
            max_rel_diff = fabs((C_dense[i*N+j]-spamm_recursive_get(i, j, C))/C_dense[i*N+j]);
          }
        }
      }
    }
  }

  printf("max diff = %e, rel. diff = %e, A(%u,%u) = %e, A_reference(%u,%u) = %e\n",
      max_diff,
      (C_dense[max_i*N+max_j] != 0.0 ? max_diff/C_dense[max_i*N+max_j] : 0.0),
      max_i, max_j, spamm_recursive_get(max_i, max_j, C),
      max_i, max_j, C_dense[max_i*N+max_j]);

  if (max_rel_diff > TEST_TOLERANCE)
  {
    printf("test failed\n");
    result = -1;
  }

  return result;
}
