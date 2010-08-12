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

#define ALIGNMENT 64
#define N_BLOCK 4

//#undef HAVE_POSIX_MEMALIGN

//#define EXTERNAL_BLAS
//#define STREAM_KERNEL_1
//#define STREAM_KERNEL_2
#define STREAM_KERNEL_3

//#define DEBUG_LOOP

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
};

struct multiply_stream_index_t
{
  unsigned long long index;

  unsigned int A_index;
  unsigned int B_index;
  unsigned int C_index;
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

/* Swap two elements of type struct multiply_stream_index_t.
 */
void
swap_multiply_stream_index (struct multiply_stream_index_t *a, struct multiply_stream_index_t *b)
{
  unsigned long long temp_index;
  unsigned int temp_block_index;

  temp_index = a->index; a->index = b->index; b->index = temp_index;

  temp_block_index = a->A_index; a->A_index = b->A_index; b->A_index = temp_block_index;
  temp_block_index = a->B_index; a->B_index = b->B_index; b->B_index = temp_block_index;
  temp_block_index = a->C_index; a->C_index = b->C_index; b->C_index = temp_block_index;
}

void
quicksort_multiply_stream_index (const unsigned int left, const unsigned int right, struct multiply_stream_index_t *stream)
{
  unsigned int i, pivot, new_pivot;
  unsigned long long index;

  if (right > left)
  {
    pivot = left+(right-left)/2;
    index = stream[pivot].index;
    swap_multiply_stream_index(&stream[pivot], &stream[right]); /* Move pivot to end. */
    new_pivot = left;
    for (i = left; i < right; i++)
    {
      if (stream[i].index <= index)
      {
        if (i != new_pivot)
        {
          swap_multiply_stream_index(&stream[i], &stream[new_pivot]);
        }
        new_pivot++;
      }
    }
    swap_multiply_stream_index(&stream[new_pivot], &stream[right]); /* Move pivot to final place. */

    if (new_pivot > left)
    {
      quicksort_multiply_stream_index(left, new_pivot-1, stream);
    }

    if (new_pivot < right)
    {
      quicksort_multiply_stream_index(new_pivot+1, right, stream);
    }
  }
}

/* Sort stream index. */
void
sort_multiply_stream_index (const unsigned int number_elements, struct multiply_stream_index_t *stream)
{
  quicksort_multiply_stream_index(0, number_elements-1, stream);
}

/* Bit-Interleaves 3 integer indices.
 */
unsigned int
interleave_3_index (const unsigned int i, const unsigned int j, const unsigned int k)
{
  unsigned int result = 0;
  unsigned int bitmask = 1;
  unsigned int addmask = 1;
  unsigned int index;

  for (index = 0; index < sizeof(unsigned int)*8; index++)
  {
    if (i & bitmask) { result |= addmask; }
    addmask <<= 1;
    if (j & bitmask) { result |= addmask; }
    addmask <<= 1;
    if (k & bitmask) { result |= addmask; }
    addmask <<= 1;
    bitmask <<= 1;
  }

  return result;
}

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

  printf("linear: ");
  for (i = 0; i < N*N; i++)
  {
    printf(" %f", A[i]);
  }
  printf("\n");
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
stream_multiply (const unsigned long long number_stream_elements,
    const float alpha,
    struct multiply_stream_t *multiply_stream)
{
  unsigned long long stream_index;
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
  A = multiply_stream[0].A_block;
  B = multiply_stream[0].B_block;
  C = multiply_stream[0].C_block;

  for (stream_index = 0; stream_index < number_stream_elements; stream_index++) {
    //A = multiply_stream[stream_index].A_block;
    //B = multiply_stream[stream_index].B_block;
    //C = multiply_stream[stream_index].C_block;

    for (i = 0; i < N_BLOCK; i++) {
      for (j = 0; j < N_BLOCK; j++) {
        for (k = 0; k < N_BLOCK; k++)
        {
          C[get_block_index(N_BLOCK, i, j)] += alpha*A[get_block_index(N_BLOCK, i, k)]*B[get_block_index(N_BLOCK, k, j)];
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
  unsigned long long stream_index, number_stream_elements;

  float alpha = 1.2;
  float beta = 0.5;

  int verify = 0;
  double max_diff;

  int print = 0;

  int random_elements = 1;

  int sort_stream = 0;

  unsigned int loop;
  unsigned long long loops = 1;
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
  struct multiply_stream_index_t *multiply_stream_index;

  unsigned int memory_access_length, old_A_index, old_B_index, old_C_index;

  int parse;
  int longindex;
  char *short_options = "hN:l:vprs";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "loops", required_argument, NULL, 'l' },
    { "verify", no_argument, NULL, 'v' },
    { "print", no_argument, NULL, 'p' },
    { "no-random", no_argument, NULL, 'r' },
    { "sort", no_argument, NULL, 's' },
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
        printf("--print       Print matrices\n");
        printf("--no-random   Full matrices with index values as opposed to random\n");
        printf("--sort        Sort stream\n");
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

      case 'p':
        print = 1;
        break;

      case 'r':
        random_elements = 0;
        break;

      case 's':
        sort_stream = 1;
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
  if (posix_memalign((void**) &A, ALIGNMENT, sizeof(float)*N*N) != 0)
  {
    printf("error allocating A\n");
    exit(1);
  }
  if (posix_memalign((void**) &B, ALIGNMENT, sizeof(float)*N*N) != 0)
  {
    printf("error allocating B\n");
    exit(1);
  }
  if (posix_memalign((void**) &C, ALIGNMENT, sizeof(float)*N*N) != 0)
  {
    printf("error allocating C\n");
    exit(1);
  }
  if (posix_memalign((void**) &D, ALIGNMENT, sizeof(float)*N*N) != 0)
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

  printf("Allocated matrices at: A = %p, B = %p, C = %p, D = %p\n", A, B, C, D);

  /* Fill matrices with random data. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      if (random_elements)
      {
        A[get_index(N, i, j)] = rand()/(float) RAND_MAX;
        B[get_index(N, i, j)] = rand()/(float) RAND_MAX;
        C[get_index(N, i, j)] = rand()/(float) RAND_MAX;
      }

      else
      {
        A[get_index(N, i, j)] = get_index(N, i, j)+1;
        B[get_index(N, i, j)] = get_index(N, i, j)+1;
        C[get_index(N, i, j)] = get_index(N, i, j)+1;
      }

      D[get_index(N, i, j)] = C[get_index(N, i, j)];
    }
  }

  if (print)
  {
    printf("A =\n");
    print_matrix(N, A);
    printf("B =\n");
    print_matrix(N, B);
    printf("C =\n");
    print_matrix(N, C);
  }

  /* Allocate stream. */
  number_stream_elements = N/N_BLOCK*N/N_BLOCK*N/N_BLOCK;
  printf("allocating stream with %llu elements\n", number_stream_elements);
  printf("1 stream element has %lu bytes\n", sizeof(struct multiply_stream_t));
  if (sizeof(struct multiply_stream_t)*number_stream_elements < 1024)
  {
    printf("multiply stream has %llu bytes\n", sizeof(struct multiply_stream_t)*number_stream_elements);
  }

  else if (sizeof(struct multiply_stream_t)*number_stream_elements < 1024*1024)
  {
    printf("multiply stream has %1.2f kB\n", sizeof(struct multiply_stream_t)*number_stream_elements/1024.);
  }

  else if (sizeof(struct multiply_stream_t)*number_stream_elements < 1024*1024*1024)
  {
    printf("multiply stream has %1.2f MB\n", sizeof(struct multiply_stream_t)*number_stream_elements/1024./1024.);
  }

  else
  {
    printf("multiply stream has %1.2f GB\n", sizeof(struct multiply_stream_t)*number_stream_elements/1024./1024./1024.);
  }

#ifdef HAVE_POSIX_MEMALIGN
  if (posix_memalign((void**) &multiply_stream, ALIGNMENT, sizeof(struct multiply_stream_t)*number_stream_elements) != 0)
  {
    printf("error allocating multiply stream\n");
    exit(1);
  }
  if (posix_memalign((void**) &multiply_stream_index, ALIGNMENT, sizeof(struct multiply_stream_index_t)*number_stream_elements) != 0)
  {
    printf("error allocating multiply stream index\n");
    exit(1);
  }
#else
  if ((multiply_stream = (struct multiply_stream_t*) malloc(sizeof(struct multiply_stream_t)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream\n");
    exit(1);
  }
  if ((multiply_stream = (struct multiply_stream_index_t*) malloc(sizeof(struct multiply_stream_index_t)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream index\n");
    exit(1);
  }
#endif

  printf("allocated multiply stream at %p\n", multiply_stream);
  printf("allocated multiply stream index at %p\n", multiply_stream_index);

  stream_index = 0;
  for (i = 0; i < N; i += N_BLOCK) {
    for (j = 0; j < N; j += N_BLOCK) {
      for (k = 0; k < N; k += N_BLOCK)
      {
        multiply_stream[stream_index].A_block = &A[get_blocked_index(N, i, k)];
        multiply_stream[stream_index].B_block = &B[get_blocked_index(N, k, j)];
        multiply_stream[stream_index].C_block = &C[get_blocked_index(N, i, j)];

        multiply_stream_index[stream_index].A_index = get_blocked_index(N, i, k)/N_BLOCK/N_BLOCK;
        multiply_stream_index[stream_index].B_index = get_blocked_index(N, k, j)/N_BLOCK/N_BLOCK;
        multiply_stream_index[stream_index].C_index = get_blocked_index(N, i, j)/N_BLOCK/N_BLOCK;
        multiply_stream_index[stream_index].index = interleave_3_index(multiply_stream_index[stream_index].A_index,
            multiply_stream_index[stream_index].B_index,
            multiply_stream_index[stream_index].C_index);
        stream_index++;
      }
    }
  }
  printf("created %llu multiply stream elements\n", stream_index);

  if (print)
  {
    printf("multiply_stream_index = ");
    for (i = 0; i < N/N_BLOCK*N/N_BLOCK*N/N_BLOCK; i++)
    {
      printf(" %llu[%u:%u:%u]", multiply_stream_index[i].index,
          multiply_stream_index[i].A_index,
          multiply_stream_index[i].B_index,
          multiply_stream_index[i].C_index);
    }
    printf("\n");
  }

  if (sort_stream)
  {
    printf("sorting stream\n");
    sort_multiply_stream_index(number_stream_elements, multiply_stream_index);

    if (print)
    {
      printf("multiply_stream_index = ");
      for (i = 0; i < N/N_BLOCK*N/N_BLOCK*N/N_BLOCK; i++)
      {
        printf(" %llu[%u:%u:%u]", multiply_stream_index[i].index,
            multiply_stream_index[i].A_index,
            multiply_stream_index[i].B_index,
            multiply_stream_index[i].C_index);
      }
      printf("\n");
    }
  }

  /* Get some statistics. */
  memory_access_length = 0;
  old_A_index = multiply_stream_index[0].A_index;
  old_B_index = multiply_stream_index[0].B_index;
  old_C_index = multiply_stream_index[0].C_index;
  for (i = 0; i < number_stream_elements; i++)
  {
    memory_access_length += abs(multiply_stream_index[i].A_index - old_A_index);
    memory_access_length += abs(multiply_stream_index[i].B_index - old_B_index);
    memory_access_length += abs(multiply_stream_index[i].C_index - old_C_index);

    old_A_index = multiply_stream_index[i].A_index;
    old_B_index = multiply_stream_index[i].B_index;
    old_C_index = multiply_stream_index[i].C_index;
  }
  printf("memory access length = %u\n", memory_access_length);

  /* Apply beta to C. */
  for (i = 0; i < N*N; i++)
  {
    C[i] *= beta;
    if (verify) { D[i] *= beta; }
  }
  printf("applied beta\n");

  printf("looping over multiply (loops = %llu)\n", loops);

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);

#ifdef HAVE_PAPI
  papi_values = (long long*) malloc(sizeof(long long)*6);
  PAPI_set_granularity(PAPI_GRN_MIN);
  PAPI_create_eventset(&papi_events);
  PAPI_add_event(papi_events, PAPI_TOT_INS);
  PAPI_add_event(papi_events, PAPI_TOT_CYC);
  PAPI_add_event(papi_events, PAPI_RES_STL);
  PAPI_add_event(papi_events, PAPI_L1_TCM);
  PAPI_add_event(papi_events, PAPI_L2_TCM);
  PAPI_add_event(papi_events, PAPI_TLB_DM);
  PAPI_start(papi_events);
#endif

#ifdef DEBUG_LOOP
  printf("debugging loop\n");
  stream_multiply(loops, alpha, multiply_stream);

#else
  for (loop = 0; loop < loops; loop++)
  {
#if defined(EXTERNAL_BLAS)
    sgemm_("N", "N", &N, &N, &N, &alpha, A, &N, B, &N, &beta, C, &N);

#else
    stream_multiply(number_stream_elements, alpha, multiply_stream);

#endif
  }
#endif

#ifdef HAVE_PAPI
  PAPI_stop(papi_events, papi_values);
  printf("[PAPI] %lli total instructions\n",    papi_values[0]);
  printf("[PAPI] %lli total cycles\n",          papi_values[1]);
  printf("[PAPI] %lli cycles stalled\n",        papi_values[2]);
  printf("[PAPI] %lli total L1 misses\n",       papi_values[3]);
  printf("[PAPI] %lli total L2 misses\n",       papi_values[4]);
  printf("[PAPI] %lli total data TLB misses\n", papi_values[5]);
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

  if (print)
  {
    printf("C =\n");
    print_matrix(N, C);
  }

  if (verify)
  {
    printf("verify\n");
    stream_multiply_verify(N, alpha, A, B, D);

    if (print)
    {
      printf("D =\n");
      print_matrix(N, D);
    }

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
