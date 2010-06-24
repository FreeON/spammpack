#include <spamm.h>
#include <math.h>
#include <stdlib.h>

int
main ()
{
  struct spamm_t A;
  struct spamm_t B;
  struct spamm_t C;

  float_t *A_dense;
  float_t *B_dense;
  float_t *C_dense;

  int i, j, k;

  double tolerance = 1e-5;

  float_t alpha = 1.5;
  float_t beta = 0.5;

  unsigned int L = 4;
  unsigned int M = 4;
  unsigned int N = 4;

  unsigned int L_block = 1;
  unsigned int M_block = 1;
  unsigned int N_block = 1;

  unsigned int M_child = 2;
  unsigned int N_child = 2;

  double fill = 1.0;

  /* Allocate memory. */
  spamm_new(L, M, L_block, M_block, M_child, N_child, 0.0, &A);
  spamm_new(M, N, M_block, N_block, M_child, N_child, 0.0, &B);
  spamm_new(L, N, L_block, N_block, M_child, N_child, 0.0, &C);

  A_dense = (float_t*) malloc(sizeof(float_t)*L*M);
  B_dense = (float_t*) malloc(sizeof(float_t)*M*N);
  C_dense = (float_t*) malloc(sizeof(float_t)*L*N);

  /* Fill matrices with random data. */
  for (i = 0; i < L; ++i) {
    for (j = 0; j < M; ++j)
    {
      A_dense[spamm_dense_index(i, j, L, M)] = (rand()/(double) RAND_MAX > (1-fill) ? rand()/(double) RAND_MAX : 0.0);
    }
  }

  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      B_dense[spamm_dense_index(i, j, M, N)] = (rand()/(double) RAND_MAX > (1-fill) ? rand()/(double) RAND_MAX : 0.0);
    }
  }

  for (i = 0; i < L; ++i) {
    for (j = 0; j < N; ++j)
    {
      C_dense[spamm_dense_index(i, j, L, N)] = (rand()/(double) RAND_MAX > (1-fill) ? rand()/(double) RAND_MAX : 0.0);
    }
  }

  /* Print some stuff for debugging. */
  //printf("A:\n");
  //spamm_print_dense(L, M, A_dense);
  //printf("B:\n");
  //spamm_print_dense(M, N, B_dense);
  //printf("C:\n");
  //spamm_print_dense(L, N, C_dense);

  /* Convert matrices. */
  spamm_dense_to_spamm(L, M, L_block, M_block, M_child, N_child, 0.0, A_dense, &A);
  spamm_dense_to_spamm(M, N, M_block, N_block, M_child, N_child, 0.0, B_dense, &B);
  spamm_dense_to_spamm(L, N, L_block, N_block, M_child, N_child, 0.0, C_dense, &C);

  /* Multiply. */
  for (i = 0; i < L; ++i) {
    for (j = 0; j < N; ++j) {
      C_dense[spamm_dense_index(i, j, L, N)] *= beta;
      for (k = 0; k < M; ++k)
      {
        C_dense[spamm_dense_index(i, j, L, N)] += alpha*A_dense[spamm_dense_index(i, k, L, M)]*B_dense[spamm_dense_index(k, j, M, N)];
      }
    }
  }

  spamm_multiply(tree, alpha, &A, &B, beta, &C);

  /* Compare. */
  for (i = 0; i < L; ++i) {
    for (j = 0; j < N; ++j)
    {
      if (fabs(C_dense[spamm_dense_index(i, j, L, N)]-spamm_get(i, j, &C)) > tolerance)
      {
        LOG_FATAL("mismatch: C_dense(%u,%u) != C(%u,%u) (%f != %f)\n", i, j, i, j, C_dense[spamm_dense_index(i, j, L, N)], spamm_get(i, j, &C));
        return -1;
      }
    }
  }

  return 0;
}
