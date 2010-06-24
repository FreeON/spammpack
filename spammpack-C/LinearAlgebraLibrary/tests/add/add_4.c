#include <spamm.h>
#include <math.h>
#include <stdlib.h>

int
main ()
{
  struct spamm_t A;
  struct spamm_t B;

  float_t *A_dense;
  float_t *B_dense;

  int i, j;

  double tolerance = 1e-5;

  float_t alpha = 1.2;
  float_t beta = 0.5;

  unsigned int M = 4;
  unsigned int N = 4;

  unsigned int M_block = 1;
  unsigned int N_block = 1;

  unsigned int M_child = 2;
  unsigned int N_child = 2;

  double fill = 0.2;

  /* Allocate memory. */
  spamm_new(M, N, M_block, N_block, M_child, N_child, 0.0, &A);
  spamm_new(M, N, M_block, N_block, M_child, N_child, 0.0, &B);

  A_dense = (float_t*) malloc(sizeof(float_t)*M*N);
  B_dense = (float_t*) malloc(sizeof(float_t)*M*N);

  /* Fill matrices with random data. */
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      A_dense[spamm_dense_index(i, j, M, N)] = (rand()/(double) RAND_MAX > (1-fill) ? rand()/(double) RAND_MAX : 0.0);
      B_dense[spamm_dense_index(i, j, M, N)] = (rand()/(double) RAND_MAX > (1-fill) ? rand()/(double) RAND_MAX : 0.0);
    }
  }

#ifdef TEST_DEBUG
  printf("A:\n");
  spamm_print_dense(M, N, A_dense);
  printf("B:\n");
  spamm_print_dense(M, N, B_dense);
#endif

  /* Convert matrices. */
  spamm_dense_to_spamm(M, N, M_block, N_block, M_child, N_child, 0.0, A_dense, &A);
  spamm_dense_to_spamm(M, N, M_block, N_block, M_child, N_child, 0.0, B_dense, &B);

#ifdef TEST_DEBUG
  printf("A (SpAMM):\n");
  spamm_print_spamm(&A);
  printf("B (SpAMM):\n");
  spamm_print_spamm(&B);
#endif

  /* Add. */
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j) {
      B_dense[spamm_dense_index(i, j, M, N)] = alpha*A_dense[spamm_dense_index(i, j, M, N)]+beta*B_dense[spamm_dense_index(i, j, M, N)];
    }
  }

#ifdef TEST_DEBUG
  printf("B:\n");
  spamm_print_dense(M, N, B_dense);
#endif

  spamm_add(alpha, &A, beta, &B);

#ifdef TEST_DEBUG
  printf("B (SpAMM):\n");
  spamm_print_spamm(&B);
#endif

  /* Compare. */
  for (i = 0; i < M; ++i) {
    for (j = 0; j < N; ++j)
    {
      if (fabs(B_dense[spamm_dense_index(i, j, M, N)]-spamm_get(i, j, &B)) > tolerance)
      {
        LOG_FATAL("mismatch: B_dense(%u,%u) != B(%u,%u) (%f != %f)\n", i, j, i, j, B_dense[spamm_dense_index(i, j, M, N)], spamm_get(i, j, &B));
        return -1;
      }
    }
  }

  return 0;
}
