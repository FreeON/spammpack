#include <spamm.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>
//#include <fenv.h>
#include <sys/time.h>
#include <sys/resource.h>

#define THRESHOLD 1e-14

//#pragma STDC FENV_ACCESS ON

enum matrix_t
{
  dense, diagonal, column_row
};

const int number_matrix_types = 3;
const char *matrix_type_name[] = { "dense", "diagonal", "column_row" };

const int number_mul_types = 3;
const char *mul_type_name[] = { "tree", "cache", "cache_redundant" };

int
main (int argc, char **argv)
{
  int parse;
  int result = 0;

  char *version_string;

  enum matrix_t matrix_type = dense;
  enum spamm_multiply_algorithm_t mul_type = cache;

  unsigned int M = 1000;
  unsigned int N = 1000;
  unsigned int K = 1000;
  unsigned int M_block = 2;
  unsigned int N_block = 2;
  unsigned int K_block = 2;

  unsigned int number_nonzeros;

  int random_elements = 1;
  int print_matrix = 0;
  int use_kahan = 0;
  int verify_spamm = 0;
  int use_spamm = 1;

  int linear_tier = -1;
  unsigned int chunksize = 100;

  //int fpe_raised;

  floating_point_t threshold = 0.0;
  floating_point_t gamma = 0.46;

  floating_point_t alpha = 1.2;
  floating_point_t beta = 0.5;

  unsigned int i, j;

#if ! defined(DGEMM)
  unsigned int k;
#endif

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime_blas, walltime_spamm;
  double usertime_blas, usertime_spamm;
  double systime_blas, systime_spamm;

  floating_point_t *A_dense;
  floating_point_t *B_dense;
  floating_point_t *C_dense;

  struct spamm_t A, B, C;

  struct spamm_tree_stats_t tree_stats;

  floating_point_t max_diff;
  unsigned int max_diff_i = 0;
  unsigned int max_diff_j = 0;

  floating_point_t temp;

  /* For the Kahan summation algorithm. */
  floating_point_t sum, error_compensation, error_compensation_max, Kahan_y, Kahan_t;

  char *short_options = "hM:N:K:1:2:3:g:e:a:b:t:m:l:c:pkvsn";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "M", required_argument, NULL, 'M' },
    { "N", required_argument, NULL, 'N' },
    { "K", required_argument, NULL, 'K' },
    { "M_block", required_argument, NULL, '1' },
    { "N_block", required_argument, NULL, '2' },
    { "K_block", required_argument, NULL, '3' },
    { "gamma", required_argument, NULL, 'g' },
    { "alpha", required_argument, NULL, 'a' },
    { "beta", required_argument, NULL, 'b' },
    { "thresh", required_argument, NULL, 'e' },
    { "type", required_argument, NULL, 't' },
    { "multype", required_argument, NULL, 'm' },
    { "linear", required_argument, NULL, 'l' },
    { "chunk", required_argument, NULL, 'c' },
    { "print", no_argument, NULL, 'p' },
    { "kahan", no_argument, NULL, 'k' },
    { "verify", no_argument, NULL, 'v' },
    { "no-spamm", no_argument, NULL, 's' },
    { "no-random", no_argument, NULL, 'n' },
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
        printf("-h                        This help\n");
        printf("-M M                      Set the number or rows of matrix A and C\n");
        printf("-N N                      Set the number of columns of matrix B and C\n");
        printf("-K K                      Set the number of columns of matrix A and the\n");
        printf("                          number or rows of matrix B\n");
        printf("--M_block M               Set the number of rows of the matrix blocks of A and C\n");
        printf("--N_block N               Set the number of columns of the matrix blocks of B and C\n");
        printf("--K_block K               Set the number of columns of the matrix blocks of A\n");
        printf("                          and the number of rows of the matrix blocks of B\n");
        printf("{ -g | --gamma } gamma    Set the decay constant gamma for exp(-gamma |i-j|)\n");
        printf("{ -a | --alpha } alpha    Set alpha in C = alpha*A*B + beta*C\n");
        printf("{ -b | --beta } beta      Set beta in C = alpha*A*B + beta*C\n");
        printf("{ -e | --thresh } eps     Set the threshold, i.e. no random number will be used below\n");
        printf("                          this threshold\n");
        printf("{ -t | --type } type      The calculation type (dense, diagonal, column_row)\n");
        printf("{ -m | --multype } type   The multiplication algorithm (tree, cache, cache_redundant)\n");
        printf("{ -l | --linear } n       Convert to linear trees at tier n\n");
        printf("{ -c | --chunk } s        Use chunks of s bytes for linear tree\n");
        printf("--print                   Print the matrices\n");
        printf("--kahan                   Use Kahan summation algorithm\n");
        printf("--verify                  Verify SpAMM by multiplying dense matrices\n");
        printf("--no-spamm                Run without any spamm operations\n");
        printf("--no-random               Fill matrices with linear index as opposed to random values\n");
        return result;
        break;

      case 'M':
        M = strtol(optarg, NULL, 10);
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'K':
        K = strtol(optarg, NULL, 10);
        break;

      case '1':
        M_block = strtol(optarg, NULL, 10);
        break;

      case '2':
        N_block = strtol(optarg, NULL, 10);
        break;

      case '3':
        K_block = strtol(optarg, NULL, 10);
        break;

      case 'g':
        gamma = strtod(optarg, NULL);
        break;

      case 'a':
        alpha = strtod(optarg, NULL);
        break;

      case 'b':
        beta = strtod(optarg, NULL);
        break;

      case 'e':
        threshold = strtod(optarg, NULL);
        break;

      case 't':
        for (i = 0; i < number_matrix_types; ++i)
        {
          if (strcasecmp(matrix_type_name[i], optarg) == 0)
          {
            matrix_type = i;
          }
        }
        break;

      case 'm':
        for (i = 0; i < number_mul_types; i++)
        {
          if (strcasecmp(mul_type_name[i], optarg) == 0)
          {
            mul_type = i;
          }
        }
        break;

      case 'l':
        linear_tier = strtol(optarg, NULL, 10);
        break;

      case 'c':
        chunksize = strtol(optarg, NULL, 10);
        break;

      case 'p':
        print_matrix = 1;
        break;

      case 'k':
        use_kahan = 1;
        break;

      case 'v':
        verify_spamm = 1;
        break;

      case 's':
        use_spamm = 0;
        break;

      case 'n':
        random_elements = 0;
        break;

      default:
        LOG2_FATAL("unknown command line option\n");
        return -1;
        break;
    }
  }

  version_string = spamm_version();
  LOG_INFO("using SpAMM version: %s\n", version_string);
  free(version_string);

  LOG_INFO("matrix mode: %s\n", matrix_type_name[matrix_type]);
  LOG_INFO("alpha = %f, beta = %f\n", alpha, beta);

  if (linear_tier >= 0)
  {
    LOG_INFO("linear tier at %i, chunks of %u bytes\n", linear_tier, chunksize);
  }

  if (use_kahan)
  {
    LOG2_INFO("using Kahan summation algorithm\n");
  }

  /* Clear floating point exceptions. */
  //feclearexcept(FE_ALL_EXCEPT);

  switch (matrix_type)
  {
    case dense:
    case column_row:
      LOG_INFO("generating random %ix%i matrix A with %ix%i blocks\n", M, K, M_block, K_block);
      LOG_INFO("generating random %ix%i matrix B with %ix%i blocks\n", K, N, K_block, N_block);
      LOG_INFO("generating random %ix%i matrix C with %ix%i blocks\n", M, N, M_block, N_block);
      break;

    case diagonal:
      if (M == N && N == K && M == K)
      {
        LOG_INFO("generating random %ix%i matrix A with %ix%i blocks, gamma = %f\n", M, K, M_block, K_block, gamma);
        LOG_INFO("generating random %ix%i matrix B with %ix%i blocks, gamma = %f\n", K, N, K_block, N_block, gamma);
        LOG_INFO("generating random %ix%i matrix C with %ix%i blocks, gamma = %f\n", M, N, M_block, N_block, gamma);
      }

      else
      {
        LOG2_FATAL("matrix dimensions do not match\n");
        exit(1);
      }
      break;

    default:
      LOG2_FATAL("unknown type\n");
      exit(1);
      break;
  }

  A_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*M*K);
  B_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*K*N);
  C_dense = (floating_point_t*) malloc(sizeof(floating_point_t)*M*N);

  /* Initialize matrices. */
  switch (matrix_type)
  {
    case dense:
      for (i = 0; i < M; i++) {
        for (j = 0; j < K; j++)
        {
          A_dense[spamm_dense_index(i, j, M, K)] = (random_elements ? rand()/(floating_point_t) RAND_MAX + threshold : spamm_dense_index(i, j, M, K));
        }
      }
      for (i = 0; i < K; i++) {
        for (j = 0; j < N; j++)
        {
          B_dense[spamm_dense_index(i, j, K, N)] = (random_elements ? rand()/(floating_point_t) RAND_MAX + threshold : spamm_dense_index(i, j, K, N));
        }
      }
      for (i = 0; i < M; i++) {
        for (j = 0; j < N; j++)
        {
          C_dense[spamm_dense_index(i, j, M, N)] = (random_elements ? rand()/(floating_point_t) RAND_MAX + threshold : spamm_dense_index(i, j, M, N));
        }
      }
      break;

    case diagonal:
      for (i = 0; i < N; i++)
      {
        A_dense[spamm_dense_index(i, i, N, N)] = rand()/(floating_point_t) RAND_MAX + threshold;
        B_dense[spamm_dense_index(i, i, N, N)] = rand()/(floating_point_t) RAND_MAX + threshold;
        C_dense[spamm_dense_index(i, i, N, N)] = rand()/(floating_point_t) RAND_MAX + threshold;
      }

      for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
        {
          if (i != j)
          {
            A_dense[spamm_dense_index(i, j, N, N)] = exp(-gamma*fabs((int) i - (int) j))*A_dense[spamm_dense_index(i, i, N, N)]*A_dense[spamm_dense_index(j, j, N, N)];
            B_dense[spamm_dense_index(i, j, N, N)] = exp(-gamma*fabs((int) i - (int) j))*B_dense[spamm_dense_index(i, i, N, N)]*B_dense[spamm_dense_index(j, j, N, N)];
            C_dense[spamm_dense_index(i, j, N, N)] = exp(-gamma*fabs((int) i - (int) j))*C_dense[spamm_dense_index(i, i, N, N)]*C_dense[spamm_dense_index(j, j, N, N)];
          }
        }
      }
      break;

    case column_row:
      for (i = 0; i < M; ++i) {
        for (j = 0; j < K; ++j)
        {
          if (j == 0)
          {
            A_dense[spamm_dense_index(i, j, M, K)] = rand()/(floating_point_t) RAND_MAX + threshold;
          }

          else
          {
            A_dense[spamm_dense_index(i, j, M, K)] = 0.0;
          }
        }
      }

      for (i = 0; i < K; ++i) {
        for (j = 0; j < N; ++j)
        {
          if (i == 0)
          {
            B_dense[spamm_dense_index(i, j, K, N)] = rand()/(floating_point_t) RAND_MAX + threshold;
          }

          else
          {
            B_dense[spamm_dense_index(i, j, K, N)] = 0.0;
          }
        }
      }

      for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j)
        {
          if (i == 0)
          {
            C_dense[spamm_dense_index(i, j, M, N)] = rand()/(floating_point_t) RAND_MAX + threshold;
          }

          else
          {
            C_dense[spamm_dense_index(i, j, M, N)] = 0.0;
          }
        }
      }
      break;
  }

  if (print_matrix)
  {
    printf("A (dense):\n");
    spamm_print_dense(M, K, A_dense);
    printf("B (dense):\n");
    spamm_print_dense(K, N, B_dense);
    printf("C (dense):\n");
    spamm_print_dense(M, N, C_dense);
  }

  number_nonzeros = 0;
  for (i = 0; i < M; i++) {
    for (j = 0; j < K; j++)
    {
      if (A_dense[spamm_dense_index(i, j, M, K)] != 0.0) { number_nonzeros++; }
    }
  }
  LOG_INFO("A has %u nonzeros, avg. of %1.2f nonzeros per row, %1.2f%% sparsity:\n",
      number_nonzeros, number_nonzeros/(double) M, (1-number_nonzeros/(double) (M*K))*100);

  number_nonzeros = 0;
  for (i = 0; i < K; i++) {
    for (j = 0; j < N; j++)
    {
      if (B_dense[spamm_dense_index(i, j, K, N)] != 0.0) { number_nonzeros++; }
    }
  }
  LOG_INFO("B has %u nonzeros, avg. of %1.2f nonzeros per row, %1.2f%% sparsity:\n",
      number_nonzeros, number_nonzeros/(double) K, (1-number_nonzeros/(double) (K*N))*100);

  number_nonzeros = 0;
  for (i = 0; i < M; i++) {
    for (j = 0; j < N; j++)
    {
      if (C_dense[spamm_dense_index(i, j, M, N)] != 0.0) { number_nonzeros++; }
    }
  }
  LOG_INFO("C has %u nonzeros, avg. of %1.2f nonzeros per row, %1.2f%% sparsity:\n",
      number_nonzeros, number_nonzeros/(double) M, (1-number_nonzeros/(double) (M*N))*100);

  if (use_spamm)
  {
    /* Convert to SpAMM. */
    spamm_new(M, K, M_block, K_block, 2, 2, 0.0, &A);
    spamm_new(K, N, K_block, N_block, 2, 2, 0.0, &B);
    spamm_new(M, N, M_block, N_block, 2, 2, 0.0, &C);

    LOG2_INFO("converting dense to spamm... ");
    gettimeofday(&start, NULL);
    spamm_dense_to_spamm(M, K, M_block, K_block, 2, 2, 0.0, A_dense, &A);
    spamm_dense_to_spamm(K, N, K_block, N_block, 2, 2, 0.0, B_dense, &B);
    spamm_dense_to_spamm(M, N, M_block, N_block, 2, 2, 0.0, C_dense, &C);
    gettimeofday(&stop, NULL);
    printf("walltime %f s\n", (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6);

    if (linear_tier >= 0)
    {
      spamm_tree_pack(linear_tier, chunksize, i_mask, &A);
      spamm_tree_pack(linear_tier, chunksize, i_mask, &B);
      spamm_tree_pack(linear_tier, chunksize, i_mask, &C);
    }

    spamm_tree_stats(&tree_stats, &A);
    LOG_INFO("A: %ix%i, padded %ix%i, depth = %u, avg. sparsity = %1.1f%%, dense blocks = %u\n", A.M, A.N, A.M_padded, A.N_padded, A.tree_depth, tree_stats.average_sparsity*100, tree_stats.number_dense_blocks);
    spamm_tree_stats(&tree_stats, &B);
    LOG_INFO("B: %ix%i, padded %ix%i, depth = %u, avg. sparsity = %1.1f%%, dense blocks = %u\n", B.M, B.N, B.M_padded, B.N_padded, B.tree_depth, tree_stats.average_sparsity*100, tree_stats.number_dense_blocks);
    spamm_tree_stats(&tree_stats, &C);
    LOG_INFO("C: %ix%i, padded %ix%i, depth = %u, avg. sparsity = %1.1f%%, dense blocks = %u\n", C.M, C.N, C.M_padded, C.N_padded, C.tree_depth, tree_stats.average_sparsity*100, tree_stats.number_dense_blocks);

    if (print_matrix)
    {
      printf("A:\n");
      spamm_print_spamm(&A);
      printf("B:\n");
      spamm_print_spamm(&B);
      printf("C:\n");
      spamm_print_spamm(&C);
    }
  }

  if (verify_spamm)
  {
#if defined(DGEMM)
    LOG2_INFO("multiplying matrix with BLAS to get reference product... ");
    gettimeofday(&start, NULL);
    getrusage(RUSAGE_SELF, &rusage_start);
    DGEMM("N", "N", &M, &N, &K, &alpha, A_dense, &M, B_dense, &K, &beta, C_dense, &M);
    getrusage(RUSAGE_SELF, &rusage_stop);
    gettimeofday(&stop, NULL);
    walltime_blas = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
    usertime_blas = (rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6;
    systime_blas = (rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6;
    printf("walltime: %f s, user time: %f s + system time: %f s = %f s\n", walltime_blas, usertime_blas, systime_blas, usertime_blas+systime_blas);
#else
    /* Clear floating point exceptions. */
    //fpe_raised = fetestexcept(FE_ALL_EXCEPT);
    //LOG_INFO("fp exceptions raised: %u = ", fpe_raised);
    //if (FE_DIVBYZERO | fpe_raised) { printf("DIVBYZERO "); }
    //if (FE_INEXACT | fpe_raised) { printf("INEXACT "); }
    //if (FE_INVALID | fpe_raised) { printf("INVALID "); }
    //if (FE_OVERFLOW | fpe_raised) { printf("OVERFLOW "); }
    //if (FE_UNDERFLOW | fpe_raised) { printf("UNDERFLOW "); }
    //printf("\n");
    //feclearexcept(FE_ALL_EXCEPT);

    LOG2_INFO("multiplying matrix directly to get reference product... ");
    gettimeofday(&start, NULL);
    getrusage(RUSAGE_SELF, &rusage_start);
    error_compensation_max = 0;
    for (i = 0; i < M; i++) {
      for (j = 0; j < N; j++) {
        C_dense[spamm_dense_index(i, j, M, N)] *= beta;

        if (use_kahan)
        {
          sum = C_dense[spamm_dense_index(i, j, M, N)];
          error_compensation = 0;
        }

        for (k = 0; k < K; k++)
        {
          if (use_kahan)
          {
            C_dense[spamm_dense_index(i, j, M, N)] += rand()*rand();
            Kahan_y = alpha*A_dense[spamm_dense_index(i, k, M, K)]*B_dense[spamm_dense_index(k, j, K, N)] - error_compensation;
            Kahan_t = sum + Kahan_y;
            error_compensation = (Kahan_t - sum) - Kahan_y;
            if (fabs(error_compensation) > fabs(error_compensation_max)) { error_compensation_max = error_compensation; }
            sum = Kahan_t;

            /* Check for floating point exceptions. */
            //feclearexcept(FE_ALL_EXCEPT);
          }

          else
          {
            C_dense[spamm_dense_index(i, j, M, N)] += alpha*A_dense[spamm_dense_index(i, k, M, K)]*B_dense[spamm_dense_index(k, j, K, N)];
          }
        }

        if (use_kahan)
        {
          C_dense[spamm_dense_index(i, j, M, N)] = sum;
        }
      }
    }
    getrusage(RUSAGE_SELF, &rusage_stop);
    gettimeofday(&stop, NULL);
    walltime_blas = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
    usertime_blas = (rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6;
    systime_blas = (rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6;
    printf("walltime: %f s, user time: %f s + system time: %f s = %f s\n", walltime_blas, usertime_blas, systime_blas, usertime_blas+systime_blas);
    if (use_kahan)
    {
      LOG_INFO("max error_compensation = %f\n", error_compensation_max);
    }
#endif

    if (print_matrix)
    {
      printf("result C (dense):\n");
      spamm_print_dense(M, N, C_dense);
    }
  }

  if (use_spamm)
  {
    LOG2_INFO("multiplying matrix with spamm\n");
    gettimeofday(&start, NULL);
    getrusage(RUSAGE_SELF, &rusage_start);
    spamm_multiply(mul_type, alpha, &A, &B, beta, &C);
    getrusage(RUSAGE_SELF, &rusage_stop);
    gettimeofday(&stop, NULL);
    walltime_spamm = (stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6;
    usertime_spamm = (rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6;
    systime_spamm = (rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6;
    LOG_INFO("product has %u nonzero elements\n", spamm_number_nonzero(&C));
    LOG_INFO("total spamm time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s\n", walltime_spamm, usertime_spamm, systime_spamm, usertime_spamm+systime_spamm);

    if (print_matrix)
    {
      printf("result C:\n");
      spamm_print_spamm(&C);
    }
  }

  if (verify_spamm && use_spamm)
  {
    /* Print out some exciting results. */
    LOG_INFO("ratio between SpAMM and BLAS: SpAMM = %f x BLAS\n", walltime_spamm/walltime_blas);

    LOG2_INFO("comparing matrices\n");
    max_diff = 0;
    for (i = 0; i < C.M; i++) {
      for (j = 0; j < C.N; j++)
      {
        if (fabs(spamm_get(i, j, &C)-C_dense[spamm_dense_index(i, j, C.M, C.N)]) > max_diff)
        {
          max_diff = fabs(spamm_get(i, j, &C)-C_dense[spamm_dense_index(i, j, C.M, C.N)]);
          max_diff_i = i;
          max_diff_j = j;
        }
      }
    }

    if (max_diff > THRESHOLD)
    {
      printf("[multiply_spamm] biggest mismatch above threshold of %e: (C[%i][%i] = %e) != (C_dense[%i][%i] = %e), |diff| = %e\n",
          THRESHOLD,
          max_diff_i, max_diff_j, spamm_get(max_diff_i, max_diff_j, &C),
          max_diff_i, max_diff_j, C_dense[spamm_dense_index(max_diff_i, max_diff_j, C.M, C.N)],
          fabs(spamm_get(max_diff_i, max_diff_j, &C)-C_dense[spamm_dense_index(max_diff_i, max_diff_j, C.M, C.N)]));
      result = 1;
    }

    if (result == 0)
    {
      LOG2_INFO("comparison ok, matrices are identical\n");
    }
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
