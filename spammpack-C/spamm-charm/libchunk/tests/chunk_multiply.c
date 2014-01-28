#include "config.h"

#include "chunk.h"
#include "lapack_interface.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>

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

  double gamma = 4.0;

  enum
  {
    full, decay, diagonal
  }
  matrix_type = full;

  short print_complexity = 0;
  short print_matrix = 0;
  short verify = 1;
  short tree_only = 0;
  short test_dense_product = 0;

  double tolerance = 0;

  int c;
  const char short_options[] = "hT:ct:N:b:pvrg:d";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "type", required_argument, NULL, 'T' },
    { "complexity", no_argument, NULL, 'c' },
    { "tolerance", required_argument, NULL, 't' },
    { "N_chunk", required_argument, NULL, 'N' },
    { "N_basic", required_argument, NULL, 'b' },
    { "print", no_argument, NULL, 'p' },
    { "no-verify", no_argument, NULL, 'v' },
    { "tree-only", no_argument, NULL, 'r' },
    { "gamma", required_argument, NULL, 'g' },
    { "dense", no_argument, NULL, 'd' },
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
        printf("{ -T | --type } TYPE      The matrix type: { full, decay, diagonal }\n");
        printf("{ -c | --complexity}      Print product complexity\n");
        printf("{ -t | --tolerance } TOL  The SpAMM tolerance\n");
        printf("{ -N | --N_chunk } N      The matrix size, N_chunk\n");
        printf("{ -b | --N_basic } N      The basic sub-matrix size, N_basic\n");
        printf("{ -p | --print }          Print matrices\n");
        printf("{ -v | --no-verify }      Do not verify result\n");
        printf("{ -r | --tree-only }      Skip basic block products\n");
        printf("{ -g | --gamma } GAMMA    The decay constant\n");
        printf("{ -d | --dense }          Multiply the matrix with a dense method for timing\n");
        exit(0);
        break;

      case 'T':
        if(strcasecmp(optarg, "full") == 0)
        {
          matrix_type = full;
        }

        else if(strcasecmp(optarg, "decay") == 0)
        {
          matrix_type = decay;
        }

        else if(strcasecmp(optarg, "diagonal") == 0)
        {
          matrix_type = diagonal;
        }

        else
        {
          printf("unknown matrix type\n");
          exit(-1);
        }
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

      case 'g':
        gamma = strtod(optarg, NULL);
        break;

      case 'd':
        test_dense_product = 1;
        break;

      default:
        printf("illegal argument\n");
        exit(-1);
        break;
    }
  }

  void *A = chunk_alloc(N_chunk, N_basic, N, 0, 0);
  void *C = chunk_alloc(N_chunk, N_basic, N, 0, 0);

  /* Fill in column-major order. */
  double *A_dense = calloc(N_chunk*N_chunk, sizeof(double));

  printf("allocated A_dense, sizeof(A_dense) = %ld bytes\n", N_chunk*N_chunk*sizeof(double));

  switch(matrix_type)
  {
    case full:
      for(int i = 0; i < N_chunk*N_chunk; i++)
      {
        A_dense[i] = rand()/(double) RAND_MAX;
      }
      break;

    case decay:
      printf("decay constant = %e\n", gamma);
      for(int i = 0; i < N_chunk; i++)
      {
        A_dense[COLUMN_MAJOR(i, i, N_chunk)] = rand()/(double) RAND_MAX;
        for(int j = 0; j < N_chunk; j++)
        {
          if(i != j)
          {
            A_dense[COLUMN_MAJOR(i, j, N_chunk)] = A_dense[COLUMN_MAJOR(i, i, N_chunk)]*exp(-fabs(i-j)/gamma);
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
  printf("random matrix, norm_2 = %e, norm_2 of chunk = %e\n",
      norm_2, chunk_get_norm(A));

  struct timespec start_time;
  clock_gettime(CLOCKTYPE, &start_time);
  chunk_multiply(tolerance, A, A, C, tree_only);
  struct timespec end_time;
  clock_gettime(CLOCKTYPE, &end_time);

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

  if(test_dense_product)
  {
    struct timespec start_time;
    clock_gettime(CLOCKTYPE, &start_time);
    double *C_dense = calloc(SQUARE(N_chunk), sizeof(double));
    double alpha = 1.0;
    double beta = 1.0;
    DGEMM("N", "N", &N_chunk, &N_chunk, &N_chunk, &alpha, A_dense, &N_chunk,
        A_dense, &N_chunk, &beta, C_dense, &N_chunk);
    struct timespec end_time;
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

  if(verify)
  {
    printf("verifying...\n");

    double *C_dense = chunk_to_dense(C);

    if(print_complexity)
    {
      int complexity = 0;
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
                norm_A += SQUARE(A_dense[(i+i_basic*N_basic)*N_chunk+(k+j_basic*N_basic)]);
                norm_B += SQUARE(A_dense[(k+i_basic*N_basic)*N_chunk+(j+j_basic*N_basic)]);
              }
            }

            if(norm_A*norm_B > tolerance)
            {
              complexity++;
            }
          }
        }
      }

      printf("product complexity = %d out of %d\n",
          complexity, CUBED(N_chunk/N_basic));
    }

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
  free(A);
  free(C);

  return 0;
}
