#include <spamm.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#define THRESHOLD 1e-14

int
main (int argc, char **argv)
{
  int parse;
  int result = 0;

  unsigned int N = 1000;
  unsigned int M_block = 100;
  unsigned int N_block = 100;

#ifdef DGEMM
  float_t alpha, beta;
#endif

  unsigned int i, j;

#if ! defined(DGEMM)
  unsigned int k;
#endif

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime_blas, walltime_spamm;
  double usertime_blas, usertime_spamm;

  float_t *A_dense;
  float_t *B_dense;
  float_t *C_dense;

  struct spamm_t A, B, C;

  float_t max_diff;
  unsigned int max_diff_i = 0;
  unsigned int max_diff_j = 0;

  char *short_options = "hN:1:2:";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "N_block", required_argument, NULL, '1' },
    { "M_block", required_argument, NULL, '2' },
    { NULL, 0, NULL, 0 }
  };
  int longindex;

  spamm_set_loglevel(info);

  while ((parse = getopt_long(argc, argv, short_options, long_options, &longindex)) != -1)
  {
    switch (parse)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("-h             This help\n");
        printf("-N N           Set the size of the matrix to N\n");
        printf("--N_block N    Set the size of the number of rows of the matrix blocks to N\n");
        printf("--M_block M    Set the size of the number of columns of the matrix blocks to M\n");
        return result;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case '1':
        M_block = strtol(optarg, NULL, 10);
        break;

      case '2':
        N_block = strtol(optarg, NULL, 10);
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
      if (j == 0)
      {
        A_dense[spamm_dense_index(i, j, N, N)] = rand()/(float_t) RAND_MAX;
      }

      else
      {
        A_dense[spamm_dense_index(i, j, N, N)] = 0.0;
      }

      if (i == 0)
      {
        B_dense[spamm_dense_index(i, j, N, N)] = rand()/(float_t) RAND_MAX;
      }

      else
      {
        B_dense[spamm_dense_index(i, j, N, N)] = 0.0;
      }

      C_dense[spamm_dense_index(i, j, N, N)] = 0.0;
    }
  }

#if defined(DGEMM)
  LOG2_INFO("multiplying matrix with BLAS to get reference product... ");
  alpha = 1.0;
  beta = 0.0;
  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  DGEMM("N", "N", &N, &N, &N, &alpha, A_dense, &N, B_dense, &N, &beta, C_dense, &N);
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_blas = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
  usertime_blas = (rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6;
  printf("walltime: %f s, user time: %f s\n", walltime_blas, usertime_blas);
#else
  LOG2_INFO("multiplying matrix directly to get reference product... ");
  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  for (i = 0; i < N; ++i) {
    for (j = 0; j < N; ++j) {
      for (k = 0; k < N; ++k)
      {
        C_dense[spamm_dense_index(i, j, N, N)] += A_dense[spamm_dense_index(i, k, N, N)]*B_dense[spamm_dense_index(k, j, N, N)];
      }
    }
  }
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_blas = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
  usertime_blas = (rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6;
  printf("walltime: %f s, user time: %f s\n", walltime_blas, usertime_blas);
#endif

#ifdef BENCHMARK_DEBUG
  printf("A (dense):\n");
  spamm_print_dense(N, N, A_dense);
  printf("B (dense):\n");
  spamm_print_dense(N, N, B_dense);
  printf("C (dense):\n");
  spamm_print_dense(N, N, C_dense);
#endif

  /* Convert to SpAMM. */
  spamm_new(N, N, M_block, N_block, 2, 2, 0.0, &A);
  spamm_new(N, N, N_block, M_block, 2, 2, 0.0, &B);
  spamm_new(N, N, M_block, M_block, 2, 2, 0.0, &C);

  LOG2_INFO("converting dense to spamm... ");
  gettimeofday(&start, NULL);
  spamm_dense_to_spamm(N, N, M_block, N_block, 2, 2, 0.0, A_dense, &A);
  spamm_dense_to_spamm(N, N, N_block, M_block, 2, 2, 0.0, B_dense, &B);
  gettimeofday(&stop, NULL);
  printf("walltime %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

#ifdef BENCHMARK_DEBUG
  printf("A:\n");
  spamm_print_spamm(&A);
  printf("B:\n");
  spamm_print_spamm(&B);
#endif

  LOG2_INFO("multiplying matrix with spamm\n");
  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  spamm_multiply(cache, 1.0, &A, &B, 1.0, &C);
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_spamm = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
  usertime_spamm = (rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6;
  LOG_INFO("product has %u nonzero elements\n", spamm_number_nonzero(&C));
  LOG_INFO("total spamm time elapsed, walltime: %f s, usertime: %f s\n", walltime_spamm, usertime_spamm);

  /* Print out some exciting results. */
  LOG_INFO("ratio between SpAMM and BLAS: SpAMM = %f x BLAS\n", walltime_spamm/walltime_blas);

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
