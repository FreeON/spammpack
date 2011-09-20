#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "spamm.h"
#include "spamm_naive.h"

#define TOLERANCE 5e-7

#ifndef SGEMM
#define SGEMM sgemm_
#endif

//#define PRINT_DEBUG

int
main ()
{
  const int N =  10;
  const int blocksize = 3;
  const float alpha = 1.2;
  const float beta = 0.5;

  float tolerance = 1e-4;

  struct spamm_naive_t *A = NULL;
  struct spamm_naive_t *B = NULL;
  struct spamm_naive_t *C = NULL;
  float *A_dense = NULL;
  float *B_dense = NULL;
  float *C_dense = NULL;

  int i, j, k;

  A = spamm_naive_new(N, N, blocksize);
  B = spamm_naive_new(N, N, blocksize);
  C = spamm_naive_new(N, N, blocksize);

  A_dense = calloc(N*N, sizeof(float));
  B_dense = calloc(N*N, sizeof(float));
  C_dense = calloc(N*N, sizeof(float));

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      A_dense[spamm_index_row_major(i, j, N, N)] = rand()/(float) RAND_MAX;
      B_dense[spamm_index_row_major(i, j, N, N)] = rand()/(float) RAND_MAX;
      C_dense[spamm_index_row_major(i, j, N, N)] = rand()/(float) RAND_MAX;

      spamm_naive_set(i, j, A_dense[spamm_index_row_major(i, j, N, N)], A);
      spamm_naive_set(i, j, B_dense[spamm_index_row_major(i, j, N, N)], B);
      spamm_naive_set(i, j, C_dense[spamm_index_row_major(i, j, N, N)], C);
    }
  }


  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (A_dense[spamm_index_row_major(i, j, N, N)] != spamm_naive_get(i, j, A))
      {
        printf("[%s:%i] element mismatch\n", __FILE__, __LINE__);
        return -1;
      }

      if (B_dense[spamm_index_row_major(i, j, N, N)] != spamm_naive_get(i, j, B))
      {
        printf("[%s:%i] element mismatch\n", __FILE__, __LINE__);
        return -1;
      }

      if (C_dense[spamm_index_row_major(i, j, N, N)] != spamm_naive_get(i, j, C))
      {
        printf("[%s:%i] element mismatch\n", __FILE__, __LINE__);
        return -1;
      }
    }
  }

#ifdef PRINT_DEBUG
  /* Print matrices for debugging. */
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

  printf("C =\n");
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
  printf("C =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", C_dense[spamm_index_row_major(i, j, N, N)]);
    }
    printf("\n");
  }
#endif

  spamm_naive_multiply(tolerance, alpha, A, B, beta, C, NULL, spamm_naive_sgemm);

#ifdef PRING_DEBUG
  printf("C(SpAMM) =\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" % f", spamm_naive_get(i, j, C));
    }
    printf("\n");
  }
#endif

  /* Compare results. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (C_dense[spamm_index_row_major(i, j, N, N)] != 0.0)
      {
        if (fabs((C_dense[spamm_index_row_major(i, j, N, N)]-spamm_naive_get(i, j, C))/C_dense[spamm_index_row_major(i, j, N, N)]) > TOLERANCE)
        {
          printf("[%s:%i] result mismatch: |(C_dense(%i,%i)-C(%i,%i))/C_dense(%i,%i)| = %e\n", __FILE__, __LINE__,
              i, j, i, j, i, j,
              fabs((C_dense[spamm_index_row_major(i, j, N, N)]-spamm_naive_get(i, j, C))/C_dense[spamm_index_row_major(i, j, N, N)]));
          return -1;
        }
      }
    }
  }

  return 0;
}
