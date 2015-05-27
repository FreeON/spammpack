#include "config.h"

#include "bcsr.h"
#include "chunk.h"
#include "lapack_interface.h"

#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <sys/mman.h>
#include <time.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQUARE(x) ((x)*(x))
#define CUBED(x) ((x)*(x)*(x))

#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))
#define COLUMN_MAJOR(i, j, N) ((i)+(j)*(N))
#define INDEX(i, j, N) COLUMN_MAJOR(i, j, N)

void
print_dense (const double *A, const int N, const char *label, ...)
{
  va_list ap;
  char expanded_label[2000];

  va_start(ap, label);
  vsnprintf(expanded_label, 2000, label, ap);
  va_end(ap);

  printf("%s\n", expanded_label);
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      printf(" % 1.3f", A[COLUMN_MAJOR(i, j, N)]);
    }
    printf("\n");
  }
}

int
main (int argc, char **argv)
{
  int N = 1024;
  int N_chunk = 1024;
  short N_chunk_set = 0;
  int N_basic = 4;
  int repeat = 1;

  /* Decay constant for large matrices. */
  double lambda = 0.995;

  enum
  {
    full, exponential_decay, algebraic_decay, diagonal, BCSR
  }
  matrix_type = full;

  char *BCSR_filename = NULL;

  short print_complexity = 0;
  short print_matrix = 0;
  short verify = 1;
  short tree_only = 0;
  short test_dense_product = 0;

  double tolerance = 0;

  int c;
  const char short_options[] = "hT:f:Cc:t:N:b:pvrl:dR:";
  const struct option long_options[] = {
    { "help",       no_argument,       NULL, 'h' },
    { "type",       required_argument, NULL, 'T' },
    { "bcsr",       required_argument, NULL, 'f' },
    { "complexity", no_argument,       NULL, 'C' },
    { "tolerance",  required_argument, NULL, 't' },
    { "N",          required_argument, NULL, 'N' },
    { "N_chunk",    required_argument, NULL, 'c' },
    { "N_basic",    required_argument, NULL, 'b' },
    { "print",      no_argument,       NULL, 'p' },
    { "no-verify",  no_argument,       NULL, 'v' },
    { "tree-only",  no_argument,       NULL, 'r' },
    { "lambda",     required_argument, NULL, 'l' },
    { "dense",      no_argument,       NULL, 'd' },
    { "repeat",     required_argument, NULL, 'R' },
    { NULL, 0, NULL, 0 }
  };

  while((c = getopt_long(argc, argv,
          short_options, long_options, NULL)) != -1)
  {
    switch(c)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("{ -h | --help }           This help\n");
        printf("{ -T | --type } TYPE      The matrix type: { full, exp_decay, "
               "alg_decay, diagonal, BCSR }\n");
        printf("{ -f | --bcsr } FILE      Load BCSR matrix from FILE\n");
        printf("{ -C | --complexity}      Print product complexity\n");
        printf("{ -t | --tolerance } TOL  The SpAMM tolerance\n");
        printf("{ -N | --N } N            The matrix size N\n");
        printf("{ -c | --N_chunk } N      The chunk size, N_chunk\n");
        printf("{ -b | --N_basic } N      The basic sub-matrix size, N_basic\n");
        printf("{ -p | --print }          Print matrices\n");
        printf("{ -v | --no-verify }      Do not verify result\n");
        printf("{ -r | --tree-only }      Skip basic block products\n");
        printf("{ -l | --lambda } LAMBDA  The decay constant\n");
        printf("{ -d | --dense }          Multiply the matrix with a dense method for timing\n");
        printf("{ -R | --repeat } N       Repeat multiply N times\n");
        exit(0);
        break;

      case 'T':
        if(strcasecmp(optarg, "full") == 0)
        {
          matrix_type = full;
        }

        else if(strcasecmp(optarg, "exp_decay") == 0)
        {
          matrix_type = exponential_decay;
        }

        else if(strcasecmp(optarg, "alg_decay") == 0)
        {
          matrix_type = algebraic_decay;
        }

        else if(strcasecmp(optarg, "diagonal") == 0)
        {
          matrix_type = diagonal;
        }

        else if(strcasecmp(optarg, "BCSR") == 0)
        {
          matrix_type = BCSR;
        }

        else
        {
          printf("unknown matrix type\n");
          exit(-1);
        }
        break;

      case 'f':
        BCSR_filename = strdup(optarg);
        break;

      case 'C':
        print_complexity = 1;
        break;

      case 't':
        tolerance = strtod(optarg, NULL);
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'c':
        N_chunk = strtol(optarg, NULL, 10);
        N_chunk_set = 1;
        break;

      case 'b':
        N_basic = strtol(optarg, NULL, 10);
        break;

      case 'p':
        print_matrix = 1;
        break;

      case 'v':
        verify = 0;
        break;

      case 'r':
        tree_only = 1;
        break;

      case 'l':
        lambda = strtod(optarg, NULL);
        break;

      case 'd':
        test_dense_product = 1;
        break;

      case 'R':
        repeat = strtol(optarg, NULL, 10);
        if(repeat < 1)
        {
          printf("repeat can not be less than 1\n");
          exit(1);
        }
        break;

      default:
        printf("illegal argument\n");
        exit(-1);
        break;
    }
  }

  printf("running: chunk_multiply");
  for(int i = 1; i < argc; i++)
  {
    printf(" %s", argv[i]);
  }
  printf("\n");

  struct bcsr_t *A_BCSR = NULL;
  double *A_dense = NULL;

  int N_matrix = 1;

  if(matrix_type == BCSR)
  {
    /* Load matrix from file before allocating chunks (we need to know the
     * size).
     */
    int tier;
    A_BCSR = bcsr_load(BCSR_filename);
    N = bcsr_get_N(A_BCSR);

    if(!N_chunk_set)
    {
      printf("storing matrix in one chunk\n");
      N_chunk = N;

      /* Zero pad the chunks. */
      for(tier = 0; ; tier++)
      {
        int N_tier = N_basic*(1 << tier);
        int remainder = N_chunk%N_tier;
        if(remainder != 0)
        {
          N_chunk += N_tier-remainder;
        }

        if(N_chunk <= N_tier)
        {
          break;
        }
      }
      printf("matrix chunk has depth %d, N_chunk = %d\n", tier, N_chunk);
    }

    N_matrix = N/N_chunk;
    if(N%N_chunk != 0) N_matrix++;

    printf("storing %dx%d chunk matrix\n", N_matrix, N_matrix);

    int M_BCSR;
    int N_BCSR;
    A_dense = bcsr_to_dense(&M_BCSR, &N_BCSR, A_BCSR);
  }

  else
  {
    /* Fill in column-major order. */
    A_dense = calloc(N*N, sizeof(double));
    N_matrix = N/N_chunk;
    if(N%N_chunk != 0) N_matrix++;

    printf("N_matrix = %d\n", N_matrix);
    printf("allocated A_dense, sizeof(A_dense) = %ld bytes\n", N*N*sizeof(double));

    switch(matrix_type)
    {
      case full:
        for(int i = 0; i < N*N; i++)
        {
          A_dense[i] = rand()/(double) RAND_MAX;
        }
        break;

      case exponential_decay:
        printf("exponential decay, lambda = %e\n", lambda);
        for(int i = 0; i < N; i++)
        {
          A_dense[COLUMN_MAJOR(i, i, N)] = 0.5+0.5*(rand()/(double) RAND_MAX);
          for(int j = i+1; j < N; j++)
          {
            A_dense[COLUMN_MAJOR(i, j, N)] = A_dense[COLUMN_MAJOR(i, i, N)]
              * exp(log(lambda)*fabs(i-j));
            A_dense[COLUMN_MAJOR(j, i, N)] = A_dense[COLUMN_MAJOR(i, j, N)];
          }
        }
        printf("A[1][N]/A[1][1] = %e\n",
               A_dense[COLUMN_MAJOR(0, N-1, N)]/A_dense[0]);
        break;

      case algebraic_decay:
        printf("algebraic decay, lambda = %e\n", lambda);
        for(int i = 0; i < N; i++)
        {
          A_dense[COLUMN_MAJOR(i, i, N)] = 0.5+0.5*(rand()/(double) RAND_MAX);
          for(int j = 0; j < N; j++)
          {
            if(i != j)
            {
              A_dense[COLUMN_MAJOR(i, j, N)] = A_dense[COLUMN_MAJOR(i, i, N)]
                / (exp(lambda*log(fabs(i-j))) + 1);
            }
          }
        }
        break;

      case diagonal:
        for(int i = 0; i < N; i++)
        {
          A_dense[COLUMN_MAJOR(i, i, N)] = rand()/(double) RAND_MAX;
        }
        break;

      default:
        printf("[FIXME] unknown matrix type\n");
        exit(-1);
        break;
    }

    if(print_matrix)
    {
      print_dense(A_dense, N, "A:");
    }
  }

  /* Zero pad matrix. */
  int N_padded = N_matrix*N_chunk;
  double *A_temp = calloc(SQUARE(N_padded), sizeof(double));
#pragma omp parallel for
  for(int i = 0; i < N; i++)
  {
    for(int j = 0; j < N; j++)
    {
      A_temp[COLUMN_MAJOR(i, j, N_padded)] = A_dense[COLUMN_MAJOR(i, j, N)];
    }
  }
  free(A_dense);
  A_dense = A_temp;

  size_t padding_count = SQUARE(N_padded)-SQUARE(N);
  printf("zero padding %dx%d BCSR matrix by %ld elements to %dx%d\n",
         N, N, padding_count, N_padded, N_padded);

  if(print_matrix)
  {
    print_dense(A_dense, N_padded, "A(padded):");
  }

  /* Convert the matrices. The chunks are stored in a N_matrix x
     N_matrix array. */
  printf("allocating matrix of %dx%d chunks\n", N_matrix, N_matrix);
  void **A = calloc(SQUARE(N_matrix), sizeof(void*));
  void **C = calloc(SQUARE(N_matrix), sizeof(void*));

  printf("allocating temporary matrix %dx%d\n", N_chunk, N_chunk);
  double *A_slice = calloc(SQUARE(N_chunk), sizeof(double));
  for(int i = 0; i < N_matrix; i++)
  {
    for(int j = 0; j < N_matrix; j++)
    {
      A[COLUMN_MAJOR(i, j, N_matrix)] = chunk_alloc(N_chunk, N_basic, N, i*N_chunk, j*N_chunk);
      C[COLUMN_MAJOR(i, j, N_matrix)] = chunk_alloc(N_chunk, N_basic, N, i*N_chunk, j*N_chunk);

      /* Set the chunk. */
      for(int i_chunk = 0; i_chunk < N_chunk; i_chunk++)
      {
        for(int j_chunk = 0; j_chunk < N_chunk; j_chunk++)
        {
          A_slice[COLUMN_MAJOR(i_chunk, j_chunk, N_chunk)] =
            A_dense[COLUMN_MAJOR(i*N_chunk+i_chunk, j*N_chunk+j_chunk, N_padded)];
        }
      }
      chunk_set(A[COLUMN_MAJOR(i, j, N_matrix)], A_slice);

      if(print_matrix)
      {
        print_dense(A_slice, N_chunk, "A_slice[%d,%d]:", i, j);

        chunk_print(A[COLUMN_MAJOR(i, j, N_matrix)], "A_slice\n");
      }
    }
  }
  free(A_slice);

  double norm_2 = 0;
  for(int i = 0; i < N_padded*N_padded; i++)
  {
    norm_2 += A_dense[i]*A_dense[i];
  }
  printf("matrix, norm_2 = %e\n", norm_2);

  if(test_dense_product)
  {
#ifdef DGEMM
    struct timespec start_time;
    struct timespec end_time;
    double *C_dense = calloc(SQUARE(N_chunk), sizeof(double));
    double alpha = 1.0;
    double beta = 1.0;
    clock_gettime(CLOCKTYPE, &start_time);
    for(int i = 0; i < repeat; i++)
    {
      DGEMM("N", "N", &N_chunk, &N_chunk, &N_chunk, &alpha, A_dense, &N_chunk,
          A_dense, &N_chunk, &beta, C_dense, &N_chunk);
    }
    clock_gettime(CLOCKTYPE, &end_time);

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        printf("done multiplying dense using %d OpenMP threads a %dx%d chunk, "
               "max. task depth %d, "
               "%1.2f seconds (%d repeats)\n",
               omp_get_num_threads(),
               N_chunk, N_chunk,
               CHUNK_TREE_MAX_TIER,
               (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
               (start_time.tv_sec+start_time.tv_nsec/1.0e9),
               repeat);
      }
    }
#else
    printf("done multiplying dense in serial a %dx%d chunk, "
           "%1.2f seconds (%d repeats)\n",
           N_chunk, N_chunk,
           (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
           (start_time.tv_sec+start_time.tv_nsec/1.0e9),
           repeat);
#endif
#else
    printf("not configured with dgemm()\n");
#endif
  }

  else
  {
    struct timespec start_time;
    struct timespec end_time;

    struct timespec individual_start_time[repeat];
    struct timespec individual_end_time[repeat];

    /* Prevent paging. */
    if(mlockall(MCL_CURRENT) != 0)
    {
      printf("can not lock memory\n");
    }

    else
    {
      printf("successfully locked memory\n");
    }

#pragma omp parallel
    {
#pragma omp single
      {
        clock_gettime(CLOCKTYPE, &start_time);
        for(int i = 0; i < repeat; i++)
        {
          clock_gettime(CLOCKTYPE, &individual_start_time[i]);
          for(int i_chunk = 0; i_chunk < N_matrix; i_chunk++)
          {
            for(int j_chunk = 0; j_chunk < N_matrix; j_chunk++)
            {
              for(int k_chunk = 0; k_chunk < N_matrix; k_chunk++)
              {
                chunk_multiply(tolerance,
                               A[COLUMN_MAJOR(i_chunk, k_chunk, N_matrix)],
                               A[COLUMN_MAJOR(k_chunk, j_chunk, N_matrix)],
                               C[COLUMN_MAJOR(i_chunk, j_chunk, N_matrix)],
                               tree_only);
                if(print_matrix)
                {
                  chunk_print(C[COLUMN_MAJOR(i_chunk, j_chunk, N_matrix)], "C:");
                }
              }
            }
          }
          clock_gettime(CLOCKTYPE, &individual_end_time[i]);
        }
        clock_gettime(CLOCKTYPE, &end_time);
      }
    }
    double total_time = ((end_time.tv_sec+end_time.tv_nsec/1.0e9)
                         -(start_time.tv_sec+start_time.tv_nsec/1.0e9))/repeat;

    double individual_duration[repeat];
    double individual_mean = 0;
    double individual_variance = 0;
    for(int i = 0; i < repeat; i++)
    {
      individual_duration[i] =
        (individual_end_time[i].tv_sec+individual_end_time[i].tv_nsec/1.0e9)-
        (individual_start_time[i].tv_sec+individual_start_time[i].tv_nsec/1.0e9);
      individual_mean += individual_duration[i];
    }
    individual_mean /= repeat;

    for(int i = 0; i < repeat; i++)
    {
      individual_variance +=
        (individual_mean-individual_duration[i])
        *(individual_mean-individual_duration[i]);
    }
    individual_variance /= repeat;

    printf("individual mean = %1.2e, variance = %1.2e\n", individual_mean, individual_variance);

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        printf("%d OpenMP threads, %dx%d matrix with %dx%d chunk and %dx%d basic blocks, "
            "tolerance %1.2e, %1.3e +- %1.3e seconds (%1.3e, %d repeats)\n",
            omp_get_num_threads(),
            N, N,
            N_chunk, N_chunk,
            N_basic, N_basic,
            tolerance,
            individual_mean,
            sqrt(individual_variance),
            total_time,
            repeat);
      }
    }
#else
    printf("serial, %dx%d matrix with %dx%d chunk and %dx%d basic blocks, "
        "tolerance %1.2e, %1.3e +- %1.3e seconds (%1.3e, %d repeats)\n",
        N, N,
        N_chunk, N_chunk,
        N_basic, N_basic,
        tolerance,
        individual_mean,
        sqrt(individual_variance),
        total_time,
        repeat);
#endif
  }

  if(print_matrix)
  {
    for(int i_chunk = 0; i_chunk < N_matrix; i_chunk++)
    {
      for(int j_chunk = 0; j_chunk < N_matrix; j_chunk++)
      {
        double *C_dense = chunk_to_dense(C[COLUMN_MAJOR(i_chunk, j_chunk, N_matrix)]);
        print_dense(C_dense, N_chunk, "C[%d,%d]:", i_chunk, j_chunk);
        free(C_dense);

        chunk_print(C[COLUMN_MAJOR(i_chunk, j_chunk, N_matrix)], "C\n");
      }
    }
  }

  if(print_complexity)
  {
    double complexity = 0;
#pragma omp parallel for reduction(+:complexity)
    for(int i = 0; i < N_padded/N_basic; i++)
    {
      for(int j = 0; j < N_padded/N_basic; j++)
      {
        for(int k = 0; k < N_padded/N_basic; k++)
        {
          double norm_A = 0;
          double norm_B = 0;

          for(int i_basic = 0; i_basic < N_basic; i_basic++)
          {
            for(int j_basic = 0; j_basic < N_basic; j_basic++)
            {
              norm_A += SQUARE(A_dense[COLUMN_MAJOR(i*N_basic+i_basic,
                                                    k*N_basic+j_basic,
                                                    N_padded)]);
              norm_B += SQUARE(A_dense[COLUMN_MAJOR(k*N_basic+i_basic,
                                                    j*N_basic+j_basic,
                                                    N_padded)]);
            }
          }

          if(norm_A*norm_B > SQUARE(tolerance))
          {
            complexity++;
          }
        }
      }
    }

    printf("product complexity = %1.0f out of %1.0f, complexity ratio = %1.3f\n",
           CUBED(N_basic)*complexity, (double) N * (double) N * (double) N,
           CUBED(N_basic)*complexity/((double) N * (double) N * (double) N));
  }

  if(verify)
  {
    printf("calculating exact product...\n");

    double *C_exact = calloc(SQUARE(N_padded), sizeof(double));

#pragma omp parallel for
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        for(int k = 0; k < N; k++)
        {
          C_exact[COLUMN_MAJOR(i, j, N_padded)] +=
            A_dense[COLUMN_MAJOR(i, k, N_padded)]*A_dense[COLUMN_MAJOR(k, j, N_padded)];
        }
      }
    }

    if(print_matrix)
    {
      print_dense(C_exact, N, "C_exact:");
    }

    double max_diff = 0;
    int max_diff_i = -1;
    int max_diff_j = -1;
    double max_diff_exact = 0;
    double max_diff_spamm = 0;

    printf("comparing exact vs. spamm...\n");
    for(int i_chunk = 0; i_chunk < N_matrix; i_chunk++)
    {
      for(int j_chunk = 0; j_chunk < N_matrix; j_chunk++)
      {
        double *C_dense = chunk_to_dense(C[COLUMN_MAJOR(i_chunk, j_chunk, N_matrix)]);

        for(int i = 0; i < N_chunk; i++)
        {
          for(int j = 0; j < N_chunk; j++)
          {
            double C_temp = C_exact[COLUMN_MAJOR(i+i_chunk*N_chunk,
                                                 j+j_chunk*N_chunk,
                                                 N_padded)];
            double diff = fabs(C_temp-C_dense[COLUMN_MAJOR(i, j, N_chunk)]);
            if(C_temp != 0)
            {
              diff = diff/C_temp;
            }

            if(diff > max_diff)
            {
              max_diff = diff;
              max_diff_i = i+i_chunk*N_chunk;
              max_diff_j = j+j_chunk*N_chunk;
              max_diff_exact = C_temp;
              max_diff_spamm = C_dense[COLUMN_MAJOR(i, j, N_chunk)];
            }
          }
        }
        free(C_dense);
      }
    }

    printf("max. mismatch at [%d][%d]: C_spamm = %e "
           " <-> C_exact = %e, rel. diff = %e\n",
           max_diff_i, max_diff_j,
           max_diff_spamm, max_diff_exact,
           max_diff);
    if(max_diff > 1e-10)
    {
      printf("max. mismatch > 1e-10\n");
      return -1;
    }
    printf("matrices are identical\n");

    free(C_exact);
  }

  free(A_dense);
  for(int i = 0; i < N_matrix; i++)
  {
    for(int j = 0; j < N_matrix; j++)
    {
      chunk_delete(&A[COLUMN_MAJOR(i, j, N_matrix)]);
      chunk_delete(&C[COLUMN_MAJOR(i, j, N_matrix)]);
    }
  }

  return 0;
}
