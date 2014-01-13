#include "config.h"

#include "chunk.h"

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <strings.h>
#include <time.h>

#define SQUARE(x) ((x)*(x))
#define CUBED(x) ((x)*(x)*(x))

int
main (int argc, char **argv)
{
  int N = 1024;
  int N_chunk = 1024;
  const int N_basic = 4;

  enum
  {
    full, diagonal
  }
  matrix_type = full;

  short print_complexity = 0;

  double tolerance = 0;

  int c;
  const char short_options[] = "hT:ct:N:";
  const struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "type", required_argument, NULL, 'T' },
    { "complexity", no_argument, NULL, 'c' },
    { "tolerance", required_argument, NULL, 't' },
    { "N_chunk", required_argument, NULL, 'N' },
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
        printf("{ -T | --type } TYPE      The matrix type: { full, diagonal }\n");
        printf("{ -c | --complexity}      Print product complexity\n");
        printf("{ -t | --tolerance } TOL  The SpAMM tolerance\n");
        printf("{ -N | --N_chunk } N      The matrix size, N_chunk\n");
        exit(0);
        break;

      case 'T':
        if(strcasecmp(optarg, "full") == 0)
        {
          matrix_type = full;
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

      default:
        printf("illegal argument\n");
        exit(-1);
        break;
    }
  }

  void *A = chunk_alloc(N_chunk, N_basic, N, 0, 0);
  void *C = chunk_alloc(N_chunk, N_basic, N, 0, 0);

  double *A_dense = calloc(N_chunk*N_chunk, sizeof(double));

  switch(matrix_type)
  {
    case full:
      for(int i = 0; i < N_chunk*N_chunk; i++)
      {
        A_dense[i] = rand()/(double) RAND_MAX;
      }
      break;

    case diagonal:
      for(int i = 0; i < N_chunk; i++)
      {
        A_dense[i*N_chunk+i] = rand()/(double) RAND_MAX;
      }
      break;

    default:
      printf("[FIXME] unknown matrix type\n");
      exit(-1);
      break;
  }

  chunk_set(A, A_dense);

  struct timespec start_time;
  clock_gettime(CLOCKTYPE, &start_time);
  chunk_multiply(tolerance, A, A, C, 0);
  struct timespec end_time;
  clock_gettime(CLOCKTYPE, &end_time);

  printf("done multiplying %dx%d chunk with %dx%d basic blocks, "
      "tolerance %1.2e, %1.2f seconds, verifying...\n",
      N_chunk, N_chunk,
      N_basic, N_basic,
      tolerance,
      (end_time.tv_sec+end_time.tv_nsec/1.0e9)-
      (start_time.tv_sec+start_time.tv_nsec/1.0e9));

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
              norm_A += SQUARE(A_dense[(i*N_basic+i_basic)*N_chunk+(k*N_basic+j_basic)]);
              norm_B += SQUARE(A_dense[(k*N_basic+i_basic)*N_chunk+(j*N_basic+j_basic)]);
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

  double C_exact;
  for(int i = 0; i < N_chunk; i++)
  {
    for(int j = 0; j < N_chunk; j++)
    {
      C_exact = 0;
      for(int k = 0; k < N_chunk; k++)
      {
        C_exact += A_dense[i*N_chunk+k]*A_dense[k*N_chunk+j];
      }

      if(fabs(C_exact-C_dense[i*N_chunk+j]) > 1e-10)
      {
        printf("mismatch C[%d][%d]\n", i, j);
        return -1;
      }
    }
  }

  printf("matrices are identical\n");

  free(A_dense);
  free(C_dense);
  free(A);
  free(C);

  return 0;
}