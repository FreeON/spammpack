#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

/* Fill matrix with random elements. */
//#define DEBUG_RANDOM_ELEMENTS

/* Print matrix for debugging. */
//#define DEBUG_PRINT_MATRIX

/* Print a lot within the multiply function. */
//#define DEBUG_PRINT

/* Randomize stream. */
//#define RANDOMIZE_STREAM

/* The blocksize. */
#define N_BLOCK 4

unsigned int
get_block_index (const unsigned i, const unsigned j)
{
  return i*N_BLOCK+j;
}

/* Convert the matrix indices i, j to a linear index for a NxN matrix. */
unsigned int
get_index (const unsigned int i, const unsigned int j, const unsigned int N)
{
  unsigned int i_block, j_block;
  unsigned int index = 0;

  i_block = i/N_BLOCK;
  j_block = j/N_BLOCK;

#ifdef RANDOMIZE_STREAM
  return N_BLOCK*(i_block*N+N_BLOCK*j_block)+get_block_index(i%N_BLOCK, j%N_BLOCK);
#else
  unsigned int bit_i, bitmask, bitset;

  bitmask = 1;
  bitset = 1;
  for (bit_i = 0; bit_i < sizeof(unsigned int)*8; bit_i++)
  {
    if (bitmask & j_block)
    {
      index |= bitset;
    }
    bitset <<= 1;

    if (bitmask & i_block)
    {
      index |= bitset;
    }
    bitset <<= 1;

    bitmask <<= 1;
  }

  return N_BLOCK*N_BLOCK*index+get_block_index(i%N_BLOCK, j%N_BLOCK);
#endif
}

/* Print a dense matrix. */
void
print_matrix (const float *A, const unsigned int N)
{
  unsigned int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf("% 1.3f ", A[get_index(i, j, N)]);
    }
    printf("\n");
  }

  printf("memory layout, linear: ");
  for (i = 0; i < N*N; ++i)
  {
    printf("% 1.3f ", A[i]);
  }
  printf("\n");
}

void
multiply_symbolic (const float alpha,
    const unsigned int N,
    unsigned int A_M_lower, unsigned int A_M_upper,
    unsigned int A_N_lower, unsigned int A_N_upper,
    const float *A,
    unsigned int B_M_lower, unsigned int B_M_upper,
    unsigned int B_N_lower, unsigned int B_N_upper,
    const float *B,
    const float beta,
    unsigned int C_M_lower, unsigned int C_M_upper,
    unsigned int C_N_lower, unsigned int C_N_upper,
    float *C,
    unsigned int *stream_index,
    unsigned int *stream_A, unsigned int *stream_B, unsigned int *stream_C)
{
  unsigned int i, j;

#ifdef DEBUG_PRINT
  printf("A[%i --> %i][%i --> %i], B[%i --> %i][%i --> %i], C[%i --> %i][%i --> %i]\n",
      A_M_lower, A_M_upper, A_N_lower, A_N_upper,
      B_M_lower, B_M_upper, B_N_lower, B_N_upper,
      C_M_lower, C_M_upper, C_N_lower, C_N_upper);
#endif

  if (A_M_upper-A_M_lower == N_BLOCK)
  {
    /* We have reached the block-level. */
#ifdef DEBUG_PRINT
    printf("multiply C[%i][%i] += % 1.3f*A[%i][%i]*B[%i][%i]",
        C_M_lower, C_N_lower,
        alpha,
        A_M_lower, A_N_lower,
        B_M_lower, B_N_lower);
    printf(" --> linear C[%i] += % 1.3f*A[%i]*B[%i]\n",
        get_index(C_M_lower, C_N_lower, N),
        alpha,
        get_index(A_M_lower, A_N_lower, N),
        get_index(B_M_lower, B_N_lower, N));
#endif
    //C[get_index(C_M_lower, C_N_lower, N)] += alpha*A[get_index(A_M_lower, A_N_lower, N)]*B[get_index(B_M_lower, B_N_lower, N)];
    stream_A[*stream_index] = get_index(A_M_lower, A_N_lower, N);
    stream_B[*stream_index] = get_index(B_M_lower, B_N_lower, N);
    stream_C[*stream_index] = get_index(C_M_lower, C_N_lower, N);
    *stream_index += 1;
  }

  else
  {
    /* Recurse down. */
#ifdef DEBUG_PRINT
    printf("1: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*0, A_M_lower+(A_M_upper-A_M_lower)/2*1,
        A_N_lower+(A_N_upper-A_N_lower)/2*0, A_N_lower+(A_N_upper-A_N_lower)/2*1,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*0, B_M_lower+(B_M_upper-B_M_lower)/2*1,
        B_N_lower+(B_N_upper-B_N_lower)/2*0, B_N_lower+(B_N_upper-B_N_lower)/2*1,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*0, C_M_lower+(C_M_upper-C_M_lower)/2*1,
        C_N_lower+(C_N_upper-C_N_lower)/2*0, C_N_lower+(C_N_upper-C_N_lower)/2*1,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("2: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*0, A_M_lower+(A_M_upper-A_M_lower)/2*1,
        A_N_lower+(A_N_upper-A_N_lower)/2*1, A_N_lower+(A_N_upper-A_N_lower)/2*2,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*1, B_M_lower+(B_M_upper-B_M_lower)/2*2,
        B_N_lower+(B_N_upper-B_N_lower)/2*0, B_N_lower+(B_N_upper-B_N_lower)/2*1,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*0, C_M_lower+(C_M_upper-C_M_lower)/2*1,
        C_N_lower+(C_N_upper-C_N_lower)/2*0, C_N_lower+(C_N_upper-C_N_lower)/2*1,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("3: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*0, A_M_lower+(A_M_upper-A_M_lower)/2*1,
        A_N_lower+(A_N_upper-A_N_lower)/2*0, A_N_lower+(A_N_upper-A_N_lower)/2*1,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*0, B_M_lower+(B_M_upper-B_M_lower)/2*1,
        B_N_lower+(B_N_upper-B_N_lower)/2*1, B_N_lower+(B_N_upper-B_N_lower)/2*2,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*0, C_M_lower+(C_M_upper-C_M_lower)/2*1,
        C_N_lower+(C_N_upper-C_N_lower)/2*1, C_N_lower+(C_N_upper-C_N_lower)/2*2,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("4: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*0, A_M_lower+(A_M_upper-A_M_lower)/2*1,
        A_N_lower+(A_N_upper-A_N_lower)/2*1, A_N_lower+(A_N_upper-A_N_lower)/2*2,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*1, B_M_lower+(B_M_upper-B_M_lower)/2*2,
        B_N_lower+(B_N_upper-B_N_lower)/2*1, B_N_lower+(B_N_upper-B_N_lower)/2*2,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*0, C_M_lower+(C_M_upper-C_M_lower)/2*1,
        C_N_lower+(C_N_upper-C_N_lower)/2*1, C_N_lower+(C_N_upper-C_N_lower)/2*2,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("5: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*1, A_M_lower+(A_M_upper-A_M_lower)/2*2,
        A_N_lower+(A_N_upper-A_N_lower)/2*0, A_N_lower+(A_N_upper-A_N_lower)/2*1,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*0, B_M_lower+(B_M_upper-B_M_lower)/2*1,
        B_N_lower+(B_N_upper-B_N_lower)/2*0, B_N_lower+(B_N_upper-B_N_lower)/2*1,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*1, C_M_lower+(C_M_upper-C_M_lower)/2*2,
        C_N_lower+(C_N_upper-C_N_lower)/2*0, C_N_lower+(C_N_upper-C_N_lower)/2*1,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("6: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*1, A_M_lower+(A_M_upper-A_M_lower)/2*2,
        A_N_lower+(A_N_upper-A_N_lower)/2*1, A_N_lower+(A_N_upper-A_N_lower)/2*2,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*1, B_M_lower+(B_M_upper-B_M_lower)/2*2,
        B_N_lower+(B_N_upper-B_N_lower)/2*0, B_N_lower+(B_N_upper-B_N_lower)/2*1,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*1, C_M_lower+(C_M_upper-C_M_lower)/2*2,
        C_N_lower+(C_N_upper-C_N_lower)/2*0, C_N_lower+(C_N_upper-C_N_lower)/2*1,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("7: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*1, A_M_lower+(A_M_upper-A_M_lower)/2*2,
        A_N_lower+(A_N_upper-A_N_lower)/2*0, A_N_lower+(A_N_upper-A_N_lower)/2*1,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*0, B_M_lower+(B_M_upper-B_M_lower)/2*1,
        B_N_lower+(B_N_upper-B_N_lower)/2*1, B_N_lower+(B_N_upper-B_N_lower)/2*2,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*1, C_M_lower+(C_M_upper-C_M_lower)/2*2,
        C_N_lower+(C_N_upper-C_N_lower)/2*1, C_N_lower+(C_N_upper-C_N_lower)/2*2,
        C,
        stream_index, stream_A, stream_B, stream_C);

#ifdef DEBUG_PRINT
    printf("8: ");
#endif
    multiply_symbolic(alpha, N,
        A_M_lower+(A_M_upper-A_M_lower)/2*1, A_M_lower+(A_M_upper-A_M_lower)/2*2,
        A_N_lower+(A_N_upper-A_N_lower)/2*1, A_N_lower+(A_N_upper-A_N_lower)/2*2,
        A,
        B_M_lower+(B_M_upper-B_M_lower)/2*1, B_M_lower+(B_M_upper-B_M_lower)/2*2,
        B_N_lower+(B_N_upper-B_N_lower)/2*1, B_N_lower+(B_N_upper-B_N_lower)/2*2,
        B,
        beta,
        C_M_lower+(C_M_upper-C_M_lower)/2*1, C_M_lower+(C_M_upper-C_M_lower)/2*2,
        C_N_lower+(C_N_upper-C_N_lower)/2*1, C_N_lower+(C_N_upper-C_N_lower)/2*2,
        C,
        stream_index, stream_A, stream_B, stream_C);
  }
}

void
multiply_stream (const float alpha, const unsigned int stream_length,
    unsigned int *stream_A, unsigned int *stream_B, unsigned int *stream_C,
    const float *A, const float *B, float *C)
{
  unsigned int i, i_block, j_block, k_block;
  float x, y, z;

  //x = 1.5;
  //y = 2.5;
  //z = 0;

  for (i = 0; i < stream_length; i++) {
    for (i_block = 0; i_block < N_BLOCK; i_block++) {
      for (j_block = 0; j_block < N_BLOCK; j_block++) {
        for (k_block = 0; k_block < N_BLOCK; k_block++)
        {
          //printf("multiplying C[%i] += alpha*A[%i]*B[%i]\n", stream_C[i], stream_A[i], stream_B[i]);
          //C[stream_C[i]+get_block_index(i_block, j_block)] += alpha*A[stream_A[i]+get_block_index(i_block, k_block)]*B[stream_B[i]+get_block_index(k_block, j_block)];
          C[stream_C[i]+i_block*N_BLOCK+j_block] += alpha*A[stream_A[i]+i_block*N_BLOCK+k_block]*B[stream_B[i]+k_block*N_BLOCK+j_block];
          //C[0] += alpha*A[0]*B[0];
          //z += alpha*x*y;
        }
      }
    }
  }

  //printf("z = %e\n", z);
}

int
main (int argc, char **argv)
{
  unsigned int N = 4;
  unsigned int N_padded;

  unsigned int temp;
  unsigned int i, j, k;

  float *A, *B, *C, *D;

  unsigned int stream_index;
  unsigned int *stream_A, *stream_B, *stream_C;

  float alpha = 1.2;
  float beta = 0.5;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;

  unsigned int repeat;
  unsigned int repeat_counter_spamm = 1;

  double walltime_spamm;
  double usertime_spamm;
  double systime_spamm;
  double flops_spamm;

  char *short_options = "hN:";
  struct option long_options[] = {
    { "N", required_argument, NULL, 'N' },
    { "repeat", required_argument, NULL, 'r' },
    { NULL, 0, NULL, 0 }
  };
  int parse;
  int longindex;

  while ((parse = getopt_long(argc, argv, short_options, long_options, &longindex)) != -1)
  {
    switch (parse)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("-h           This help\n");
        printf("-N N         Test NxN matrices\n");
        printf("--repeat N   Repeat SpAMM test N times\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'r':
        repeat_counter_spamm = strtol(optarg, NULL, 10);
        break;

      default:
        printf("unknown command line option\n");
        return -1;
        break;
    }
  }

  N_padded = (int) N_BLOCK*pow(2, (int) ceil(log(N/(double) N_BLOCK)/log(2)));
  if (N_padded != N)
  {
    printf("N needs to be a power of 2, next larger matrix with %ix%i blocks would be %ix%i\n", N_BLOCK, N_BLOCK, N_padded, N_padded);
    return -1;
  }

  printf("testing %ix%i matrices with a blocksize of %ix%i\n", N, N, N_BLOCK, N_BLOCK);

#ifdef RANDOMIZE_STREAM
  printf("randomizing stream\n");
#endif

  /* Create dense NxN matrix. */
  A = (float*) malloc(sizeof(float)*N*N);
  B = (float*) malloc(sizeof(float)*N*N);
  C = (float*) malloc(sizeof(float)*N*N);
  D = (float*) malloc(sizeof(float)*N*N);

  /* Fill with random stuff. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
#ifdef DEBUG_RANDOM_ELEMENTS
      A[get_index(i, j, N)] = rand()/(float) RAND_MAX;
      B[get_index(i, j, N)] = rand()/(float) RAND_MAX;
      C[get_index(i, j, N)] = rand()/(float) RAND_MAX;
#else
      A[get_index(i, j, N)] = get_index(i, j, N);
      B[get_index(i, j, N)] = get_index(i, j, N);
      C[get_index(i, j, N)] = get_index(i, j, N);
#endif
    }
  }

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      D[get_index(i, j, N)] = C[get_index(i, j, N)];
    }
  }

#ifdef DEBUG_PRINT_MATRIX
  printf("A =\n");
  print_matrix(A, N);
  printf("B =\n");
  print_matrix(B, N);
  printf("C =\n");
  print_matrix(C, N);
  printf("D =\n");
  print_matrix(D, N);
#endif

  /* Allocate multiply stream. */
  stream_A = (unsigned int*) malloc(sizeof(unsigned int)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK);
  stream_B = (unsigned int*) malloc(sizeof(unsigned int)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK);
  stream_C = (unsigned int*) malloc(sizeof(unsigned int)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK);

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  for (repeat = 0; repeat < repeat_counter_spamm; repeat++)
  {
    for (i = 0; i < N; i++) {
      for (j = 0; j < N; j++) {
        C[get_index(i, j, N)] *= beta;
      }
    }
    stream_index = 0;
    multiply_symbolic(alpha, N, 0, N, 0, N, A, 0, N, 0, N, B, beta, 0, N, 0, N, C, &stream_index, stream_A, stream_B, stream_C);
  }
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_spamm = ((stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6)/repeat_counter_spamm;
  usertime_spamm = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/repeat_counter_spamm;
  systime_spamm = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/repeat_counter_spamm;
  flops_spamm = ((double) N)*((double) N)*(2.0*N+1.0)/walltime_spamm;
  if (flops_spamm < 1000*1000*1000)
  {
    printf("symbolic: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Mflop/s\n",
        walltime_spamm, usertime_spamm, systime_spamm, usertime_spamm+systime_spamm, flops_spamm/1000./1000.);
  }
  else
  {
    printf("symbolic: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Gflop/s\n",
        walltime_spamm, usertime_spamm, systime_spamm, usertime_spamm+systime_spamm, flops_spamm/1000./1000./1000.);
  }

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  for (repeat = 0; repeat < repeat_counter_spamm; repeat++)
  {
    multiply_stream(alpha, N/N_BLOCK*N/N_BLOCK*N/N_BLOCK, stream_A, stream_B, stream_C, A, B, C);
  }
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_spamm = ((stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6)/repeat_counter_spamm;
  usertime_spamm = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/repeat_counter_spamm;
  systime_spamm = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/repeat_counter_spamm;
  flops_spamm = ((double) N)*((double) N)*(2.0*N+1.0)/walltime_spamm;
  if (flops_spamm < 1000*1000*1000)
  {
    printf("numeric: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Mflop/s\n",
        walltime_spamm, usertime_spamm, systime_spamm, usertime_spamm+systime_spamm, flops_spamm/1000./1000.);
  }
  else
  {
    printf("numeric: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Gflop/s\n",
        walltime_spamm, usertime_spamm, systime_spamm, usertime_spamm+systime_spamm, flops_spamm/1000./1000./1000.);
  }

  printf("multiplying by hand for verification\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      D[get_index(i, j, N)] *= beta;
      for (k = 0; k < N; k++)
      {
        D[get_index(i, j, N)] += alpha*A[get_index(i, k, N)]*B[get_index(k, j, N)];
      }
    }
  }

#ifdef DEBUG_PRINT_MATRIX
  printf("C =\n");
  print_matrix(C, N);
  printf("D =\n");
  print_matrix(D, N);
#endif

  /* Compare result. */
  printf("verifying result\n");
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (fabs(C[get_index(i, j, N)]-D[get_index(i, j, N)]) > 1e-10)
      {
        printf("mismatch C[%i][%i] != D[%i][%i] (%e != %e, |diff| = %e)\n",
            i, j, i, j, C[get_index(i, j, N)], D[get_index(i, j, N)],
            fabs(C[get_index(i, j, N)]-D[get_index(i, j, N)]));
        return -1;
      }
    }
  }

  free(A);
  free(B);
  free(C);
  free(D);
  free(stream_A);
  free(stream_B);
  free(stream_C);

  return 0;
}
