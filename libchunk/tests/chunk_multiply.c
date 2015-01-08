#include "config.h"

#include "bcsr.h"
#include "chunk.h"
#include "lapack_interface.h"

#include <getopt.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <time.h>
#include <sys/mman.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define SQUARE(x) ((x)*(x))
#define CUBED(x) ((x)*(x)*(x))

#define ROW_MAJOR(i, j, N) ((i)*(N)+(j))
#define COLUMN_MAJOR(i, j, N) ((i)+(j)*(N))
#define INDEX(i, j, N) COLUMN_MAJOR(i, j, N)

int
main (int argc, char **argv)
{
  int N = 1024;
  int N_chunk = 1024;
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
  const char short_options[] = "hT:f:ct:N:b:pvrl:dR:";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "type", required_argument, NULL, 'T' },
    { "bcsr", required_argument, NULL, 'f' },
    { "complexity", no_argument, NULL, 'c' },
    { "tolerance", required_argument, NULL, 't' },
    { "N_chunk", required_argument, NULL, 'N' },
    { "N_basic", required_argument, NULL, 'b' },
    { "print", no_argument, NULL, 'p' },
    { "no-verify", no_argument, NULL, 'v' },
    { "tree-only", no_argument, NULL, 'r' },
    { "lambda", required_argument, NULL, 'l' },
    { "dense", no_argument, NULL, 'd' },
    { "repeat", required_argument, NULL, 'R' },
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
        printf("{ -T | --type } TYPE      The matrix type: { full, exp_decay, alg_decay, diagonal, BCSR }\n");
        printf("{ -f | --bcsr } FILE      Load BCSR matrix from FILE\n");
        printf("{ -c | --complexity}      Print product complexity\n");
        printf("{ -t | --tolerance } TOL  The SpAMM tolerance\n");
        printf("{ -N | --N_chunk } N      The matrix size, N_chunk\n");
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

      case 'c':
        print_complexity = 1;
        break;

      case 't':
        tolerance = strtod(optarg, NULL);
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        N_chunk = N;
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

  if(matrix_type == BCSR)
  {
    /* Load matrix from file before allocating chunks (we need to know the
     * size).
     */
    A_BCSR = bcsr_load(BCSR_filename);
    N = bcsr_get_N(A_BCSR);

    /* Zero pad the matrix. */
    for(int tier = 0; ; tier++)
    {
      //printf("N = %d, tier = %d\n", N, tier);
      int N_tier = N_basic*(1 << tier);
      int remainder = N%N_tier;
      if(remainder != 0)
      {
        N += N_tier-remainder;
      }

      if(N <= N_tier)
      {
        break;
      }
    }
    N_chunk = N;
  }

  void *A = chunk_alloc(N_chunk, N_basic, N, 0, 0);
  void *C = chunk_alloc(N_chunk, N_basic, N, 0, 0);

  double *A_dense = NULL;

  if(matrix_type == BCSR)
  {
    int M_BCSR;
    int N_BCSR;
    A_dense = bcsr_to_dense(&M_BCSR, &N_BCSR, A_BCSR);

    /* Zero pad matrix. */
    A_dense = realloc(A_dense, sizeof(double)*N_chunk*N_chunk);
    size_t padding_count = N_chunk*N_chunk-M_BCSR*N_BCSR;
    printf("zero padding %dx%d BCSR matrix by %ld elements to N_chunk = %d\n",
        M_BCSR, N_BCSR, padding_count, N_chunk);
    memset(
        (void*) ((intptr_t) A_dense + sizeof(double)*M_BCSR*N_BCSR),
        0,
        sizeof(double)*padding_count
        );
  }

  else
  {
    /* Fill in column-major order. */
    A_dense = calloc(N_chunk*N_chunk, sizeof(double));

    printf("allocated A_dense, sizeof(A_dense) = %ld bytes\n", N_chunk*N_chunk*sizeof(double));

    switch(matrix_type)
    {
      case full:
        for(int i = 0; i < N_chunk*N_chunk; i++)
        {
          A_dense[i] = rand()/(double) RAND_MAX;
        }
        break;

      case exponential_decay:
        printf("exponential decay, lambda = %e\n", lambda);
        for(int i = 0; i < N_chunk; i++)
        {
          A_dense[COLUMN_MAJOR(i, i, N_chunk)] = 0.5+0.5*(rand()/(double) RAND_MAX);
          for(int j = i+1; j < N_chunk; j++)
          {
            A_dense[COLUMN_MAJOR(i, j, N_chunk)] = A_dense[COLUMN_MAJOR(i, i, N_chunk)] * exp(log(lambda)*fabs(i-j));
            A_dense[COLUMN_MAJOR(j, i, N_chunk)] = A_dense[COLUMN_MAJOR(i, j, N_chunk)];
          }
        }
        printf("A[1][N]/A[1][1] = %e\n", A_dense[COLUMN_MAJOR(0, N_chunk-1, N_chunk)]/A_dense[0]);
        break;

      case algebraic_decay:
        printf("algebraic decay, lambda = %e\n", lambda);
        for(int i = 0; i < N_chunk; i++)
        {
          A_dense[COLUMN_MAJOR(i, i, N_chunk)] = 0.5+0.5*(rand()/(double) RAND_MAX);
          for(int j = 0; j < N_chunk; j++)
          {
            if(i != j)
            {
              A_dense[COLUMN_MAJOR(i, j, N_chunk)] = A_dense[COLUMN_MAJOR(i, i, N_chunk)]
                / (exp(lambda*log(fabs(i-j))) + 1);
            }
          }
        }
        break;

      case diagonal:
        for(int i = 0; i < N_chunk; i++)
        {
          A_dense[COLUMN_MAJOR(i, i, N_chunk)] = rand()/(double) RAND_MAX;
        }
        break;

      default:
        printf("[FIXME] unknown matrix type\n");
        exit(-1);
        break;
    }

    if(print_matrix)
    {
      printf("A:\n");
      for(int i = 0; i < N; i++)
      {
        for(int j = 0; j < N; j++)
        {
          printf(" % 1.3f", A_dense[COLUMN_MAJOR(i, j, N_chunk)]);
        }
        printf("\n");
      }
    }
  }

  /* Set the chunk. */
  chunk_set(A, A_dense);

  if(print_matrix)
  {
    chunk_print(A, "A\n");
  }

  double norm_2 = 0;
  for(int i = 0; i < N*N; i++)
  {
    norm_2 += A_dense[i]*A_dense[i];
  }
  printf("matrix, norm_2 = %e, norm_2 of chunk = %e\n",
      norm_2, chunk_get_norm_2(A));

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
            "%1.2f seconds\n",
            omp_get_num_threads(),
            N_chunk, N_chunk,
            (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
            (start_time.tv_sec+start_time.tv_nsec/1.0e9));
      }
    }
#else
    printf("done multiplying dense in serial a %dx%d chunk, "
        "%1.2f seconds\n",
        N_chunk, N_chunk,
        (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
        (start_time.tv_sec+start_time.tv_nsec/1.0e9));
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

    clock_gettime(CLOCKTYPE, &start_time);
    for(int i = 0; i < repeat; i++)
    {
      clock_gettime(CLOCKTYPE, &individual_start_time[i]);
      chunk_multiply(tolerance, A, A, C, tree_only);
      clock_gettime(CLOCKTYPE, &individual_end_time[i]);
    }
    clock_gettime(CLOCKTYPE, &end_time);

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

    printf("individual mean = %1.2f, variance = %1.2f\n", individual_mean, individual_variance);

#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp master
      {
        printf("done multiplying using %d OpenMP threads a %dx%d chunk with %dx%d basic blocks, "
            "tolerance %1.2e, %1.2f seconds\n",
            omp_get_num_threads(),
            N_chunk, N_chunk,
            N_basic, N_basic,
            tolerance,
            (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
            (start_time.tv_sec+start_time.tv_nsec/1.0e9));
      }
    }
#else
    printf("done multiplying in serial a %dx%d chunk with %dx%d basic blocks, "
        "tolerance %1.2e, %1.2f seconds\n",
        N_chunk, N_chunk,
        N_basic, N_basic,
        tolerance,
        (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
        (start_time.tv_sec+start_time.tv_nsec/1.0e9));
#endif
  }

  if(print_matrix)
  {
    double *C_dense = chunk_to_dense(C);

    printf("C:\n");
    for(int i = 0; i < N; i++)
    {
      for(int j = 0; j < N; j++)
      {
        printf(" % 1.3f", C_dense[COLUMN_MAJOR(i, j, N_chunk)]);
      }
      printf("\n");
    }

    free(C_dense);

    chunk_print(C, "C\n");
  }

  if(print_complexity)
  {
    int complexity = 0;
#pragma omp parallel for reduction(+:complexity)
    for(int i = 0; i < N_chunk/N_basic; i++)
    {
      for(int j = 0; j < N_chunk/N_basic; j++)
      {
        for(int k = 0; k < N_chunk/N_basic; k++)
        {
          double norm_A = 0;
          double norm_B = 0;

          for(int i_basic = 0; i_basic < N_basic; i_basic++)
          {
            for(int j_basic = 0; j_basic < N_basic; j_basic++)
            {
              norm_A += SQUARE(A_dense[COLUMN_MAJOR(i*N_basic+i_basic, k*N_basic+j_basic, N_chunk)]);
              norm_B += SQUARE(A_dense[COLUMN_MAJOR(k*N_basic+i_basic, j*N_basic+j_basic, N_chunk)]);
            }
          }

          if(norm_A*norm_B > SQUARE(tolerance))
          {
            complexity++;
          }
        }
      }
    }

    printf("product complexity = %d out of %d, complexity ratio = %1.3f\n",
        complexity, CUBED(N_chunk/N_basic),
        complexity/(double) CUBED(N_chunk/N_basic));
  }

  if(verify)
  {
    printf("verifying...\n");

    double *C_dense = chunk_to_dense(C);
    double *C_exact = calloc(SQUARE(N_chunk), sizeof(double));

    for(int i = 0; i < N_chunk; i++)
    {
      for(int j = 0; j < N_chunk; j++)
      {
        for(int k = 0; k < N_chunk; k++)
        {
          C_exact[COLUMN_MAJOR(i, j, N_chunk)] +=
            A_dense[COLUMN_MAJOR(i, k, N_chunk)]*A_dense[COLUMN_MAJOR(k, j, N_chunk)];
        }
      }
    }

    if(print_matrix)
    {
      printf("C_exact:\n");
      for(int i = 0; i < N; i++)
      {
        for(int j = 0; j < N; j++)
        {
          printf(" % 1.3f", C_exact[COLUMN_MAJOR(i, j, N_chunk)]);
        }
        printf("\n");
      }
    }

    for(int i = 0; i < N_chunk; i++)
    {
      for(int j = 0; j < N_chunk; j++)
      {
        if(fabs(C_exact[COLUMN_MAJOR(i, j, N_chunk)]-C_dense[COLUMN_MAJOR(i, j, N_chunk)]) > 1e-10)
        {
          printf("mismatch C[%d][%d] = %e "
              " <-> C_exact = %e\n",
              i, j, C_dense[COLUMN_MAJOR(i, j, N_chunk)], C_exact[COLUMN_MAJOR(i, j, N_chunk)]);

          return -1;
        }
      }
    }

    printf("matrices are identical\n");

    free(C_exact);
    free(C_dense);
  }

  free(A_dense);
  chunk_delete(&A);
  chunk_delete(&C);

  return 0;
}
