#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

/* Fill matrix with random elements. */
#define DEBUG_RANDOM_ELEMENTS

/* Print matrix for debugging. */
//#define DEBUG_PRINT_MATRIX

/* Print a lot within the multiply function. */
//#define DEBUG_PRINT

/* Randomize stream. */
#define RANDOMIZE_STREAM

/* The blocksize. */
#define N_BLOCK 4

/* Use the C-code block kernel. */
//#define BLOCK_MULTIPLY_HANDCODED

/* Use the external assembly kernel. */
#define EXTERNAL_SGEMM_KERNEL_1

#if defined(EXTERNAL_SGEMM_KERNEL_1)
void
sgemm_kernel_1 (const float alpha,
    const float *restrict A,
    const float *restrict B,
    const float beta,
    float *restrict C);
#endif

/* Convert the matrix indices i, j to a linear index for a NxN matrix. */
unsigned int
get_index (const unsigned int i, const unsigned int j, const unsigned int N)
{
  unsigned int index = 0;

#ifdef RANDOMIZE_STREAM
  return i*N+j;
#else
  unsigned int bit_i, bitmask, bitset;

  bitmask = 1;
  bitset = 1;
  for (bit_i = 0; bit_i < sizeof(unsigned int)*8; bit_i++)
  {
    if (bitmask & j)
    {
      index |= bitset;
    }
    bitset <<= 1;

    if (bitmask & i)
    {
      index |= bitset;
    }
    bitset <<= 1;

    bitmask <<= 1;
  }

  return index;
#endif
}

/* Print a dense matrix. */
void
print_matrix (float (*A)[N_BLOCK][N_BLOCK], const unsigned int N)
{
  unsigned int i_block, j_block;
  unsigned int i, j;
  float *A_linear;

  for (i = 0; i < N/N_BLOCK; i++) {
    for (i_block = 0; i_block < N_BLOCK; i_block++) {
      for (j = 0; j < N/N_BLOCK; j++) {
        for (j_block = 0; j_block < N_BLOCK; j_block++)
        {
          printf("% 1.4f ", A[get_index(i, j, N/N_BLOCK)][i_block][j_block]);
        }
      }
      printf("\n");
    }
  }

  printf("memory layout, linear: ");
  A_linear = &(A[0][0][0]);
  for (i = 0; i < N*N; ++i)
  {
    printf("% 1.3f ", A_linear[i]);
  }
  printf("\n");
}

void
multiply_symbolic (const float alpha,
    const unsigned int N,
    unsigned int A_M_lower, unsigned int A_M_upper,
    unsigned int A_N_lower, unsigned int A_N_upper,
    float (*A)[N_BLOCK][N_BLOCK],
    unsigned int B_M_lower, unsigned int B_M_upper,
    unsigned int B_N_lower, unsigned int B_N_upper,
    float (*B)[N_BLOCK][N_BLOCK],
    const float beta,
    unsigned int C_M_lower, unsigned int C_M_upper,
    unsigned int C_N_lower, unsigned int C_N_upper,
    float (*C)[N_BLOCK][N_BLOCK],
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

  if (A_M_upper-A_M_lower == 1)
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
multiply_stream (const float alpha, const float beta,
    const unsigned int stream_length,
    unsigned int *stream_A, unsigned int *stream_B, unsigned int *stream_C,
    float (*A)[N_BLOCK][N_BLOCK],
    float (*B)[N_BLOCK][N_BLOCK],
    float (*C)[N_BLOCK][N_BLOCK])
{
  unsigned int i, i_block, j_block, k_block;
  unsigned int i_A, i_B, i_C;
  float x, y, z;

  for (i = 0; i < stream_length; i++)
  {
#if defined(BLOCK_MULTIPLY_HANDCODED)
    i_A = stream_A[i];
    i_B = stream_B[i];
    i_C = stream_C[i];

    C[i_C][0][0] += alpha*(
        A[i_A][0][0]*B[i_B][0][0] +
        A[i_A][0][1]*B[i_B][1][0] +
        A[i_A][0][2]*B[i_B][2][0] +
        A[i_A][0][3]*B[i_B][3][0]);

    C[i_C][0][1] += alpha*(
        A[i_A][0][0]*B[i_B][0][1] +
        A[i_A][0][1]*B[i_B][1][1] +
        A[i_A][0][2]*B[i_B][2][1] +
        A[i_A][0][3]*B[i_B][3][1]);

    C[i_C][0][2] += alpha*(
        A[i_A][0][0]*B[i_B][0][2] +
        A[i_A][0][1]*B[i_B][1][2] +
        A[i_A][0][2]*B[i_B][2][2] +
        A[i_A][0][3]*B[i_B][3][2]);

    C[i_C][0][3] += alpha*(
        A[i_A][0][0]*B[i_B][0][3] +
        A[i_A][0][1]*B[i_B][1][3] +
        A[i_A][0][2]*B[i_B][2][3] +
        A[i_A][0][3]*B[i_B][3][3]);

    C[i_C][1][0] += alpha*(
        A[i_A][1][0]*B[i_B][0][0] +
        A[i_A][1][1]*B[i_B][1][0] +
        A[i_A][1][2]*B[i_B][2][0] +
        A[i_A][1][3]*B[i_B][3][0]);

    C[i_C][1][1] += alpha*(
        A[i_A][1][0]*B[i_B][0][1] +
        A[i_A][1][1]*B[i_B][1][1] +
        A[i_A][1][2]*B[i_B][2][1] +
        A[i_A][1][3]*B[i_B][3][1]);

    C[i_C][1][2] += alpha*(
        A[i_A][1][0]*B[i_B][0][2] +
        A[i_A][1][1]*B[i_B][1][2] +
        A[i_A][1][2]*B[i_B][2][2] +
        A[i_A][1][3]*B[i_B][3][2]);

    C[i_C][1][3] += alpha*(
        A[i_A][1][0]*B[i_B][0][3] +
        A[i_A][1][1]*B[i_B][1][3] +
        A[i_A][1][2]*B[i_B][2][3] +
        A[i_A][1][3]*B[i_B][3][3]);

    C[i_C][2][0] += alpha*(
        A[i_A][2][0]*B[i_B][0][0] +
        A[i_A][2][1]*B[i_B][1][0] +
        A[i_A][2][2]*B[i_B][2][0] +
        A[i_A][2][3]*B[i_B][3][0]);

    C[i_C][2][1] += alpha*(
        A[i_A][2][0]*B[i_B][0][1] +
        A[i_A][2][1]*B[i_B][1][1] +
        A[i_A][2][2]*B[i_B][2][1] +
        A[i_A][2][3]*B[i_B][3][1]);

    C[i_C][2][2] += alpha*(
        A[i_A][2][0]*B[i_B][0][2] +
        A[i_A][2][1]*B[i_B][1][2] +
        A[i_A][2][2]*B[i_B][2][2] +
        A[i_A][2][3]*B[i_B][3][2]);

    C[i_C][2][3] += alpha*(
        A[i_A][2][0]*B[i_B][0][3] +
        A[i_A][2][1]*B[i_B][1][3] +
        A[i_A][2][2]*B[i_B][2][3] +
        A[i_A][2][3]*B[i_B][3][3]);

    C[i_C][3][0] += alpha*(
        A[i_A][3][0]*B[i_B][0][0] +
        A[i_A][3][1]*B[i_B][1][0] +
        A[i_A][3][2]*B[i_B][2][0] +
        A[i_A][3][3]*B[i_B][3][0]);

    C[i_C][3][1] += alpha*(
        A[i_A][3][0]*B[i_B][0][1] +
        A[i_A][3][1]*B[i_B][1][1] +
        A[i_A][3][2]*B[i_B][2][1] +
        A[i_A][3][3]*B[i_B][3][1]);

    C[i_C][3][2] += alpha*(
        A[i_A][3][0]*B[i_B][0][2] +
        A[i_A][3][1]*B[i_B][1][2] +
        A[i_A][3][2]*B[i_B][2][2] +
        A[i_A][3][3]*B[i_B][3][2]);

    C[i_C][3][3] += alpha*(
        A[i_A][3][0]*B[i_B][0][3] +
        A[i_A][3][1]*B[i_B][1][3] +
        A[i_A][3][2]*B[i_B][2][3] +
        A[i_A][3][3]*B[i_B][3][3]);

#elif defined(EXTERNAL_SGEMM_KERNEL_1)
    i_A = stream_A[i];
    i_B = stream_B[i];
    i_C = stream_C[i];

    //printf("A[%u] =\n", i_A);
    //print_matrix((float (*)[4][4]) A[i_A], 4);
    //printf("B[%u] =\n", i_B);
    //print_matrix((float (*)[4][4]) B[i_B], 4);
    //printf("C[%u] =\n", i_C);
    //print_matrix((float (*)[4][4]) C[i_C], 4);
    sgemm_kernel_1(alpha, (float*) &A[i_A], (float*) &B[i_B], 1.0, (float*) &C[i_C]);
    //printf("C[%u] =\n", i_C);
    //print_matrix((float (*)[4][4]) C[i_C], 4);
    //exit(0);

#endif
  }
}

int
main (int argc, char **argv)
{
  unsigned int N = 4;
  unsigned int N_padded;

  unsigned int temp;
  unsigned int i, j, k;
  unsigned int i_block, j_block, k_block;

  float (*A)[N_BLOCK][N_BLOCK];
  float (*B)[N_BLOCK][N_BLOCK];
  float (*C)[N_BLOCK][N_BLOCK];
  float (*D)[N_BLOCK][N_BLOCK];

  unsigned int stream_index;
  unsigned int *stream_A, *stream_B, *stream_C;

  float alpha = 1.2;
  float beta = 0.5;

  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;

  unsigned int repeat;
  unsigned int repeat_counter_spamm_symbolic = 1;
  unsigned int repeat_counter_apply_beta = 1;
  unsigned int repeat_counter_spamm_numeric = 1;

  double walltime_spamm_symbolic;
  double usertime_spamm_symbolic;
  double systime_spamm_symbolic;
  double flops_spamm_symbolic;

  double walltime_apply_beta;
  double usertime_apply_beta;
  double systime_apply_beta;
  double flops_apply_beta;

  double walltime_spamm_numeric;
  double usertime_spamm_numeric;
  double systime_spamm_numeric;
  double flops_spamm_numeric;

  double flops_spamm_total;

  char *short_options = "hN:";
  struct option long_options[] = {
    { "N", required_argument, NULL, 'N' },
    { "repeat-symbolic", required_argument, NULL, '1' },
    { "repeat-numeric", required_argument, NULL, '2' },
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
        printf("-h                    This help\n");
        printf("-N N                  Test NxN matrices\n");
        printf("--repeat-symbolic N   Repeat symbolic part of SpAMM test N times\n");
        printf("--repeat-numeric N    Repeat symbolic part of SpAMM test N times\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case '1':
        repeat_counter_spamm_symbolic = strtol(optarg, NULL, 10);
        break;

      case '2':
        repeat_counter_spamm_numeric = strtol(optarg, NULL, 10);
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
    N = N_padded;
  }

  printf("testing %ix%i matrices with a blocksize of %ix%i\n", N, N, N_BLOCK, N_BLOCK);

#ifdef RANDOMIZE_STREAM
  printf("randomizing stream\n");
#endif

  /* Create dense NxN matrix. */
  A = (float (*)[N_BLOCK][N_BLOCK]) malloc(sizeof(float[N_BLOCK][N_BLOCK])*N/N_BLOCK*N/N_BLOCK);
  B = (float (*)[N_BLOCK][N_BLOCK]) malloc(sizeof(float[N_BLOCK][N_BLOCK])*N/N_BLOCK*N/N_BLOCK);
  C = (float (*)[N_BLOCK][N_BLOCK]) malloc(sizeof(float[N_BLOCK][N_BLOCK])*N/N_BLOCK*N/N_BLOCK);
  D = (float (*)[N_BLOCK][N_BLOCK]) malloc(sizeof(float[N_BLOCK][N_BLOCK])*N/N_BLOCK*N/N_BLOCK);

  /* Fill with random stuff. */
  for (i = 0; i < N/N_BLOCK; i++) {
    for (i_block = 0; i_block < N_BLOCK; i_block++) {
      for (j = 0; j < N/N_BLOCK; j++) {
        for (j_block = 0; j_block < N_BLOCK; j_block++)
        {
#ifdef DEBUG_RANDOM_ELEMENTS
          A[get_index(i, j, N/N_BLOCK)][i_block][j_block] = rand()/(float) RAND_MAX;
          B[get_index(i, j, N/N_BLOCK)][i_block][j_block] = rand()/(float) RAND_MAX;
          C[get_index(i, j, N/N_BLOCK)][i_block][j_block] = rand()/(float) RAND_MAX;
#else
          A[get_index(i, j, N/N_BLOCK)][i_block][j_block] = N_BLOCK*N_BLOCK*get_index(i, j, N/N_BLOCK)+i_block*N_BLOCK+j_block;
          B[get_index(i, j, N/N_BLOCK)][i_block][j_block] = N_BLOCK*N_BLOCK*get_index(i, j, N/N_BLOCK)+i_block*N_BLOCK+j_block;
          C[get_index(i, j, N/N_BLOCK)][i_block][j_block] = N_BLOCK*N_BLOCK*get_index(i, j, N/N_BLOCK)+i_block*N_BLOCK+j_block;
#endif
        }
      }
    }
  }

  for (i = 0; i < N/N_BLOCK; i++) {
    for (i_block = 0; i_block < N_BLOCK; i_block++) {
      for (j = 0; j < N/N_BLOCK; j++) {
        for (j_block = 0; j_block < N_BLOCK; j_block++)
        {
          D[get_index(i, j, N/N_BLOCK)][i_block][j_block] = C[get_index(i, j, N/N_BLOCK)][i_block][j_block];
        }
      }
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
  for (repeat = 0; repeat < repeat_counter_spamm_symbolic; repeat++)
  {
    stream_index = 0;
    multiply_symbolic(alpha, N/N_BLOCK, 0, N/N_BLOCK, 0, N/N_BLOCK, A, 0, N/N_BLOCK, 0, N/N_BLOCK, B, beta, 0, N/N_BLOCK, 0, N/N_BLOCK, C, &stream_index, stream_A, stream_B, stream_C);
  }
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_spamm_symbolic = ((stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6)/repeat_counter_spamm_symbolic;
  usertime_spamm_symbolic = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/repeat_counter_spamm_symbolic;
  systime_spamm_symbolic = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/repeat_counter_spamm_symbolic;
  flops_spamm_symbolic = ((double) N)*((double) N)*(2.0*N+1.0)/usertime_spamm_symbolic;
  if (flops_spamm_symbolic < 1000*1000*1000)
  {
    printf("symbolic: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Mflop/s\n",
        walltime_spamm_symbolic, usertime_spamm_symbolic, systime_spamm_symbolic, usertime_spamm_symbolic+systime_spamm_symbolic, flops_spamm_symbolic/1000./1000.);
  }
  else
  {
    printf("symbolic: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Gflop/s\n",
        walltime_spamm_symbolic, usertime_spamm_symbolic, systime_spamm_symbolic, usertime_spamm_symbolic+systime_spamm_symbolic, flops_spamm_symbolic/1000./1000./1000.);
  }

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  for (repeat = 0; repeat < repeat_counter_apply_beta; repeat++)
  {
    for (i = 0; i < N/N_BLOCK; i++) {
      for (i_block = 0; i_block < N_BLOCK; i_block++) {
        for (j = 0; j < N/N_BLOCK; j++) {
          for (j_block = 0; j_block < N_BLOCK; j_block++)
          {
            C[get_index(i, j, N/N_BLOCK)][i_block][j_block] *= beta;
          }
        }
      }
    }
  }
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_apply_beta = ((stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6)/repeat_counter_apply_beta;
  usertime_apply_beta = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/repeat_counter_apply_beta;
  systime_apply_beta = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/repeat_counter_apply_beta;
  flops_apply_beta = ((double) N)*((double) N)*(2.0*N+1.0)/usertime_apply_beta;
  if (flops_apply_beta < 1000*1000*1000)
  {
    printf("apply beta: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Mflop/s\n",
        walltime_apply_beta, usertime_apply_beta, systime_apply_beta, usertime_apply_beta+systime_apply_beta, flops_apply_beta/1000./1000.);
  }
  else
  {
    printf("apply beta: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Gflop/s\n",
        walltime_apply_beta, usertime_apply_beta, systime_apply_beta, usertime_apply_beta+systime_apply_beta, flops_apply_beta/1000./1000./1000.);
  }

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);
  for (repeat = 0; repeat < repeat_counter_spamm_numeric; repeat++)
  {
    multiply_stream(alpha, beta, N/N_BLOCK*N/N_BLOCK*N/N_BLOCK, stream_A, stream_B, stream_C, A, B, C);
  }
  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime_spamm_numeric = ((stop.tv_sec-start.tv_sec)+(stop.tv_usec-start.tv_usec)/(double) 1e6)/repeat_counter_spamm_numeric;
  usertime_spamm_numeric = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/repeat_counter_spamm_numeric;
  systime_spamm_numeric = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/repeat_counter_spamm_numeric;
  flops_spamm_numeric = ((double) N)*((double) N)*(2.0*N+1.0)/usertime_spamm_numeric;
  if (flops_spamm_numeric < 1000*1000*1000)
  {
    printf("numeric: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Mflop/s\n",
        walltime_spamm_numeric, usertime_spamm_numeric, systime_spamm_numeric, usertime_spamm_numeric+systime_spamm_numeric, flops_spamm_numeric/1000./1000.);
  }
  else
  {
    printf("numeric: time elapsed, walltime: %f s, usertime: %f s + system time: %f s = %f s = %1.2f Gflop/s\n",
        walltime_spamm_numeric, usertime_spamm_numeric, systime_spamm_numeric, usertime_spamm_numeric+systime_spamm_numeric, flops_spamm_numeric/1000./1000./1000.);
  }

  flops_spamm_total = ((double) N)*((double) N)*(2.0*N+1.0)/(usertime_spamm_numeric+usertime_apply_beta+usertime_spamm_symbolic);
  if (flops_spamm_total < 1000*1000*1000)
  {
    printf("total: %1.2f Mflop/s\n", flops_spamm_total/1000./1000.);
  }
  else
  {
    printf("total: %1.2f Gflop/s\n", flops_spamm_total/1000./1000./1000.);
  }

  printf("multiplying by hand for verification\n");
  for (i = 0; i < N/N_BLOCK; i++) {
    for (i_block = 0; i_block < N_BLOCK; i_block++) {
      for (j = 0; j < N/N_BLOCK; j++) {
        for (j_block = 0; j_block < N_BLOCK; j_block++)
        {
          D[get_index(i, j, N/N_BLOCK)][i_block][j_block] *= beta;
          for (k = 0; k < N/N_BLOCK; k++) {
            for (k_block = 0; k_block < N_BLOCK; k_block++)
            {
              D[get_index(i, j, N/N_BLOCK)][i_block][j_block] += alpha*A[get_index(i, k, N/N_BLOCK)][i_block][k_block]*B[get_index(k, j, N/N_BLOCK)][k_block][j_block];
            }
          }
        }
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
  for (i = 0; i < N/N_BLOCK; i++) {
    for (i_block = 0; i_block < N_BLOCK; i_block++) {
      for (j = 0; j < N/N_BLOCK; j++) {
        for (j_block = 0; j_block < N_BLOCK; j_block++)
        {
          if (fabs(C[get_index(i, j, N/N_BLOCK)][i_block][j_block]-D[get_index(i, j, N/N_BLOCK)][i_block][j_block]) > 1e-10)
          {
            printf("mismatch C[%i][%i] != D[%i][%i] (%e != %e, |diff| = %e)\n",
                i, j, i, j, C[get_index(i, j, N/N_BLOCK)][i_block][j_block], D[get_index(i, j, N/N_BLOCK)][i_block][j_block],
                fabs(C[get_index(i, j, N/N_BLOCK)][i_block][j_block]-D[get_index(i, j, N/N_BLOCK)][i_block][j_block]));
            return -1;
          }
        }
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
