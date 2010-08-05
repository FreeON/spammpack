#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_PAPI
#include <papi.h>
#endif

#define CACHELINE_SIZE 64
#define N_BLOCK 4

//#undef HAVE_POSIX_MEMALIGN

//#define EXTERNAL_BLAS
//#define STREAM_KERNEL_1
//#define STREAM_KERNEL_2
#define STREAM_KERNEL_3

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
};

#ifdef STREAM_KERNEL_1
void
stream_kernel_1 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_2
void
stream_kernel_2 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_3
void
stream_kernel_3 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

/* Returns the Morton-ordered index of a NxN matrix. */
unsigned int
get_Morton_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return i*N+j;
}

/* Returns the linear index of a submatrix block in a (N*N_BLOCK)x(N*N_BLOCK)
 * matrix. The blocks can be arranged in any order, the stream_kernel() does
 * not need to know what order this is.
 */
unsigned int
get_block_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return get_Morton_index(N, i, j);
}

/* Returns linear index in a NxN matrix of the submatrix block. */
unsigned int
get_blocked_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  unsigned int i_block, j_block;

  i_block = (i-i%N_BLOCK)/N_BLOCK;
  j_block = (j-j%N_BLOCK)/N_BLOCK;

  return get_Morton_index(N/N_BLOCK, i_block, j_block)*N_BLOCK*N_BLOCK;
}

/* Returns the linear index of a matrix element in a blocked NxN matrix. */
unsigned int
get_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return get_blocked_index(N, i, j)+get_block_index(N_BLOCK, i%N_BLOCK, j%N_BLOCK);
}

void
print_matrix (unsigned int N, float *restrict A)
{
  unsigned int i, j;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      printf(" %f", A[get_index(N, i, j)]);
    }
    printf("\n");
  }
}

void
stream_multiply_verify (const unsigned int N, const float alpha,
    const float *restrict A, const float *restrict B, float *restrict C)
{
  unsigned int stream_index;
  unsigned int i, j, k;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++) {
      for (k = 0; k < N; k++)
      {
        C[get_index(N, i, j)] += alpha*A[get_index(N, i, k)]*B[get_index(N, k, j)];
      }
    }
  }
}

void
stream_multiply (const unsigned int number_stream_elements, const float alpha,
    struct multiply_stream_t *multiply_stream)
{
  unsigned int stream_index;
  unsigned int i, j, k;

  float *A, *B, *C;

#if defined(STREAM_KERNEL_1)
  stream_kernel_1(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_2)
  stream_kernel_2(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_3)
  stream_kernel_3(number_stream_elements, alpha, multiply_stream);

  //for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
  //{
  //  A = multiply_stream[stream_index].A_block;
  //  B = multiply_stream[stream_index].B_block;
  //  C = multiply_stream[stream_index].C_block;

  //  A[0] = 0;
  //}

#else
  /* Multiply stream. */
  for (stream_index = 0; stream_index < number_stream_elements; stream_index++) {
    for (i = 0; i < N_BLOCK; i++) {
      for (j = 0; j < N_BLOCK; j++) {
        for (k = 0; k < N_BLOCK; k++) {
          multiply_stream[stream_index].C_block[get_block_index(N_BLOCK, i, j)] +=
            alpha*multiply_stream[stream_index].A_block[get_block_index(N_BLOCK, i, k)]
            *multiply_stream[stream_index].B_block[get_block_index(N_BLOCK, k, j)];
        }
      }
    }
  }

#endif
}

double
compare_matrix (const unsigned int N, const float *restrict A, const float *restrict B)
{
  unsigned int i, j;
  double max_diff = 0;

  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (fabs(A[get_index(N, i, j)]-B[get_index(N, i, j)]) > max_diff)
      {
        max_diff = fabs(A[get_index(N, i, j)]-B[get_index(N, i, j)]);
      }
    }
  }

  return max_diff;
}

int
main (int argc, char **argv)
{
  unsigned int N = 1024;
  unsigned int N_padded = N;

  unsigned int i, j, k;
  unsigned int stream_index, number_stream_elements;

  float alpha = 1.2;
  float beta = 0.5;

  int verify = 0;
  double max_diff;

  unsigned int loop;
  unsigned int loops = 1;
  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime, flops;

#ifdef HAVE_PAPI
  int papi_events = PAPI_NULL;
  long long *papi_values;
#endif

  float *A, *B, *C, *D;
  float *A_blocked, *B_blocked, *C_blocked;

  struct multiply_stream_t *multiply_stream;

  int parse;
  int longindex;
  char *short_options = "hN:l:v";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "loops", required_argument, NULL, 'l' },
    { "verify", no_argument, NULL, 'v' },
    { NULL, 0, NULL, 0 }
  };

  /* Read command line. */
  while ((parse = getopt_long(argc, argv, short_options, long_options, &longindex)) != -1)
  {
    switch (parse)
    {
      case 'h':
        printf("Usage:\n");
        printf("\n");
        printf("-h            This help\n");
        printf("-N N          Use NxN matrix blocks\n");
        printf("--loops N     Repeat each multiply N times\n");
        printf("--verify      Verify result\n");
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loops = strtol(optarg, NULL, 10);
        break;

      case 'v':
        verify = 1;
        break;

      default:
        printf("unknown command line argument\n");
        return -1;
        break;
    }
  }

#ifdef HAVE_PAPI
  /* Do some PAPI. */
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
  {
    printf("can not initialize PAPI\n");
    exit(1);
  }
#endif

#if defined(EXTERNAL_BLAS)
  printf("external blas\n");

#elif defined(STREAM_KERNEL_1)
  printf("using stream_kernel_1\n");

#elif defined(STREAM_KERNEL_2)
  printf("using stream_kernel_2\n");

#elif defined(STREAM_KERNEL_3)
  printf("using stream_kernel_3\n");

#else
  printf("using C stream kernel\n");

#endif

  /* Check input. */
  N_padded = (int) N_BLOCK*pow(2, (int) ceil(log(N/(double) N_BLOCK)/log(2)));
  if (N_padded != N)
  {
    printf("N needs to be a power of 2, next larger matrix with %ix%i blocks would be %ix%i\n", N_BLOCK, N_BLOCK, N_padded, N_padded);
    N = N_padded;
  }

  printf("testing %ix%i matrices with a blocksize of %ix%i\n", N, N, N_BLOCK, N_BLOCK);

  /* Allocate matrices. */
#ifdef HAVE_POSIX_MEMALIGN
  if (posix_memalign((void**) &A, CACHELINE_SIZE, sizeof(float)*N*N) != 0)
  {
    printf("error allocating A\n");
    exit(1);
  }
  if (posix_memalign((void**) &B, CACHELINE_SIZE, sizeof(float)*N*N) != 0)
  {
    printf("error allocating B\n");
    exit(1);
  }
  if (posix_memalign((void**) &C, CACHELINE_SIZE, sizeof(float)*N*N) != 0)
  {
    printf("error allocating C\n");
    exit(1);
  }
  if (posix_memalign((void**) &D, CACHELINE_SIZE, sizeof(float)*N*N) != 0)
  {
    printf("error allocating D\n");
    exit(1);
  }
#else
  printf("can not allocated aligned memory for matrices\n");
  if ((A = (float*) malloc(sizeof(float)*N*N)) == NULL)
  {
    printf("error allocating A\n");
    exit(1);
  }
  if ((B = (float*) malloc(sizeof(float)*N*N)) == NULL)
  {
    printf("error allocating B\n");
    exit(1);
  }
  if ((C = (float*) malloc(sizeof(float)*N*N)) == NULL)
  {
    printf("error allocating C\n");
    exit(1);
  }
  if ((D = (float*) malloc(sizeof(float)*N*N)) == NULL)
  {
    printf("error allocating D\n");
    exit(1);
  }
#endif

  /* Fill matrices with random data. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      //A[get_index(N, i, j)] = rand()/(float) RAND_MAX;
      //B[get_index(N, i, j)] = rand()/(float) RAND_MAX;
      //C[get_index(N, i, j)] = rand()/(float) RAND_MAX;

      A[get_index(N, i, j)] = get_index(N, i, j)+1;
      B[get_index(N, i, j)] = get_index(N, i, j)+1;
      C[get_index(N, i, j)] = get_index(N, i, j)+1;

      D[get_index(N, i, j)] = C[get_index(N, i, j)];
    }
  }

  /* Allocate stream. */
  number_stream_elements = N/N_BLOCK*N/N_BLOCK*N/N_BLOCK;
  printf("allocating stream with %u elements\n", number_stream_elements);
  printf("1 stream element has %lu bytes\n", sizeof(struct multiply_stream_t));
  if (sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK < 1024)
  {
    printf("multiply stream has %lu bytes\n", sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK);
  }

  else if (sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK < 1024*1024)
  {
    printf("multiply stream has %1.2f kB\n", sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK/1024.);
  }

  else if (sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK < 1024*1024*1024)
  {
    printf("multiply stream has %1.2f MB\n", sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK/1024./1024.);
  }

  else
  {
    printf("multiply stream has %1.2f GB\n", sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK/1024./1024./1024.);
  }

#ifdef HAVE_POSIX_MEMALIGN
  if (posix_memalign((void**) &multiply_stream, CACHELINE_SIZE, sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK) != 0)
  {
    printf("error allocating multiply stream\n");
    exit(1);
  }
#else
  if ((multiply_stream = (struct multiply_stream_t*) malloc(sizeof(struct multiply_stream_t)*N/N_BLOCK*N/N_BLOCK*N/N_BLOCK)) == NULL)
  {
    printf("error allocating multiply stream\n");
    exit(1);
  }
#endif

  /* Multiply. */
  stream_index = 0;
  for (i = 0; i < N; i += N_BLOCK) {
    for (j = 0; j < N; j += N_BLOCK) {
      for (k = 0; k < N; k += N_BLOCK)
      {
        multiply_stream[stream_index].A_block = &A[get_blocked_index(N, i, k)];
        multiply_stream[stream_index].B_block = &B[get_blocked_index(N, k, j)];
        multiply_stream[stream_index].C_block = &C[get_blocked_index(N, i, j)];

        stream_index++;
      }
    }
  }
  printf("copied %u multiply stream elements\n", stream_index);

  /* Apply beta to C. */
  for (i = 0; i < N*N; i++)
  {
    C[i] *= beta;
    if (verify) { D[i] *= beta; }
  }
  printf("applied beta\n");

  printf("looping over multiply (loops = %u)\n", loops);

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);

#ifdef HAVE_PAPI
  papi_values = (long long*) malloc(sizeof(long long)*5);
  PAPI_set_granularity(PAPI_GRN_MIN);
  PAPI_create_eventset(&papi_events);
  PAPI_add_event(papi_events, PAPI_TOT_INS);
  PAPI_add_event(papi_events, PAPI_TOT_CYC);
  PAPI_add_event(papi_events, PAPI_L1_TCM);
  PAPI_add_event(papi_events, PAPI_L2_TCM);
  PAPI_add_event(papi_events, PAPI_TLB_DM);
  PAPI_start(papi_events);
#endif

  for (loop = 0; loop < loops; loop++)
  {
#if defined(EXTERNAL_BLAS)
    sgemm_("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C, &N);

#else
    stream_multiply(number_stream_elements, alpha, multiply_stream);

#endif
  }

#ifdef HAVE_PAPI
  PAPI_stop(papi_events, papi_values);
  printf("[PAPI] %lli total instructions\n",    papi_values[0]);
  printf("[PAPI] %lli total cycles\n",          papi_values[1]);
  printf("[PAPI] %lli total L1 misses\n",       papi_values[2]);
  printf("[PAPI] %lli total L2 misses\n",       papi_values[3]);
  printf("[PAPI] %lli total data TLB misses\n", papi_values[4]);
  PAPI_cleanup_eventset(papi_events);
  PAPI_destroy_eventset(&papi_events);
  free(papi_values);
#endif

  getrusage(RUSAGE_SELF, &rusage_stop);
  gettimeofday(&stop, NULL);
  walltime = (stop.tv_sec-start.tv_sec+(stop.tv_usec-start.tv_usec)/1.0e6)/loops;
  usertime = ((rusage_stop.ru_utime.tv_sec-rusage_start.ru_utime.tv_sec)+(rusage_stop.ru_utime.tv_usec-rusage_start.ru_utime.tv_usec)/(double) 1e6)/loops;
  systime = ((rusage_stop.ru_stime.tv_sec-rusage_start.ru_stime.tv_sec)+(rusage_stop.ru_stime.tv_usec-rusage_start.ru_stime.tv_usec)/(double) 1e6)/loops;
  flops = ((double) N)*((double) N)*(2.*N+1.)/usertime;
  if (flops < 1000*1000*1000)
  {
    printf("performance: total walltime %f s, usertime %f s, systime %f s, per iteration walltime %e s, usertime %e s, systime %e s, %1.2f Mflop/s\n",
        walltime*loops, usertime*loops, systime*loops,
        walltime, usertime, systime,
        flops/1000./1000.);
  }

  else
  {
    printf("performance: total walltime %f s, usertime %f s, systime %f s, per iteration walltime %e s, usertime %e s, systime %e s, %1.2f Gflop/s\n",
        walltime*loops, usertime*loops, systime*loops,
        walltime, usertime, systime,
        flops/1000./1000./1000.);
  }

  if (verify)
  {
    printf("verify\n");
    stream_multiply_verify(N, alpha, A, B, D);

    max_diff = compare_matrix(N, C, D);
    printf("max diff = %e\n", max_diff);
  }

  /* Free memory. */
  free(A);
  free(B);
  free(C);
  free(D);
  free(multiply_stream);
}
