#include <spamm.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>

#define THRESHOLD 1e-14

int
main (int argc, char **argv)
{
  unsigned int N = 1000;
  unsigned int M_block = 100;
  unsigned int N_block = 100;

  int parse;
  int result = 0;

#ifdef DGEMM
  float_t alpha, beta;
#endif

  unsigned int i, j;

#if ! defined(DGEMM)
  int k;
#endif

  struct timeval start, stop;
  double timing_blas, timing_spamm;

  float_t *A_dense;
  float_t *B_dense;
  float_t *C_dense;

  struct spamm_t A, B, C;

  float_t max_diff;
  unsigned int max_diff_i = 0;
  unsigned int max_diff_j = 0;

  spamm_set_loglevel(info);

  while ((parse = getopt(argc, argv, "hN:")) != -1)
  {
    switch (parse)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("-h      This help\n");
        printf("-N N    Set the size of the matrix to N\n");
        return result;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      default:
        LOG2_FATAL("unknown command line option\n");
        return -1;
        break;
    }
  }
  LOG_INFO("generating 2 random %ix%i matrices with %ix%i blocks\n", N, N, M_block, N_block);

  A_dense = (float_t*) malloc(sizeof(float_t)*N*N);
  B_dense = (float_t*) malloc(sizeof(float_t)*N*N);
  C_dense = (float_t*) malloc(sizeof(float_t)*N*N);


  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j)
    {
      A_dense[spamm_dense_index(i, j, N, N)] = rand()/(float_t) RAND_MAX;
      B_dense[spamm_dense_index(i, j, N, N)] = rand()/(float_t) RAND_MAX;
      C_dense[spamm_dense_index(i, j, N, N)] = 0.0;
    }
  }

#if defined(DGEMM)
  LOG2_INFO("multiplying matrix with BLAS to get reference product... ");
  alpha = 1.0;
  beta = 0.0;
  gettimeofday(&start, NULL);
  DGEMM("N", "N", &N, &N, &N, &alpha, A_dense, &N, B_dense, &N, &beta, C_dense, &N);
  gettimeofday(&stop, NULL);
  timing_blas = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
  printf("time elapsed: %f s\n", timing_blas);
#else
  LOG2_INFO("multiplying matrix directly to get reference product... ");
  gettimeofday(&start, NULL);
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k)
      {
        C_dense[spamm_dense_index(i, j, N, N)] += A_dense[spamm_dense_index(i, k, N, N)]*B_dense[spamm_dense_index(k, j, N, N)];
      }
    }
  }
  gettimeofday(&stop, NULL);
  timing_blas = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
  printf("time elapsed: %f s\n", timing_blas);
#endif

  /* Convert to SpAMM. */
  spamm_new(N, N, M_block, N_block, 2, 2, 0.0, &A);
  spamm_new(N, N, M_block, N_block, 2, 2, 0.0, &B);
  spamm_new(N, N, M_block, N_block, 2, 2, 0.0, &C);

  LOG2_INFO("converting dense to spamm... ");
  gettimeofday(&start, NULL);
  spamm_dense_to_spamm(N, N, M_block, N_block, 2, 2, 0.0, A_dense, &A);
  spamm_dense_to_spamm(N, N, M_block, N_block, 2, 2, 0.0, B_dense, &B);
  gettimeofday(&stop, NULL);
  printf("time elapsed: %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

  LOG2_INFO("multiplying matrix with spamm\n");
  gettimeofday(&start, NULL);
  spamm_multiply(cache, 1.0, &A, &B, 1.0, &C);
  gettimeofday(&stop, NULL);
  timing_spamm = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
  LOG_INFO("product has %u nonzero elements\n", spamm_number_nonzero(&C));
  LOG_INFO("total spamm time elapsed: %f s\n", timing_spamm);

  /* Print out some exciting results. */
  LOG_INFO("ratio between SpAMM and BLAS: %f\n", timing_spamm/timing_blas);

  LOG2_INFO("comparing matrices\n");
  max_diff = 0;
  for (i = 0; i < A.M; ++i) {
    for (j = 0; j < A.N; ++j)
    {
      if (fabs(spamm_get(i, j, &C)-C_dense[spamm_dense_index(i, j, A.M, A.N)]) > max_diff)
      {
        max_diff = fabs(spamm_get(i, j, &C)-C_dense[spamm_dense_index(i, j, A.M, A.N)]);
        max_diff_i = i;
        max_diff_j = j;
      }
    }
  }

  if (max_diff > THRESHOLD)
  {
    printf("[multiply_spamm] biggest mismatch above threshold of %e: (A2[%i][%i] = %e) != (A2_dense[%i][%i] = %e), |diff| = %e\n",
        THRESHOLD,
        max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &C),
        max_diff_i, max_diff_j, C_dense[spamm_dense_index(max_diff_i, max_diff_j, A.M, A.N)],
        fabs(spamm_get(max_diff_i, max_diff_j, &C)-C_dense[spamm_dense_index(max_diff_i, max_diff_j, A.M, A.N)]));
    result = 1;
  }

  /* Free memory. */
  free(A_dense);
  free(B_dense);
  free(C_dense);

  spamm_delete(&A);
  spamm_delete(&B);
  spamm_delete(&C);

  return result;
}
