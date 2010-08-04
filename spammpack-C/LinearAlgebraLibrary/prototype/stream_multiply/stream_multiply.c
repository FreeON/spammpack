#include "config.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define CACHELINE_SIZE 64
#define N_BLOCK 4

//#undef HAVE_POSIX_MEMALIGN

#define STREAM_KERNEL_1

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

unsigned int
get_Morton_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return i*N+j;
}

unsigned int
get_block_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return get_Morton_index(N, i, j);
}

unsigned int
get_blocked_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  unsigned int i_block, j_block;

  i_block = (i-i%N_BLOCK)/N_BLOCK;
  j_block = (j-j%N_BLOCK)/N_BLOCK;

  return get_Morton_index(N/N_BLOCK, i_block, j_block)*N_BLOCK*N_BLOCK;
}

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
stream_multiply (const unsigned int number_stream_elements, const float alpha, struct multiply_stream_t *multiply_stream)
{
  unsigned int stream_index;
  unsigned int i, j, k;

#if defined(STREAM_KERNEL_1)
  stream_kernel_1(number_stream_elements, alpha, multiply_stream);

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

int
main ()
{
  unsigned int N = 8;
  unsigned int N_padded = N;

  unsigned int i, j, k;
  unsigned int stream_index, number_stream_elements;

  float alpha = 1.2;
  float beta = 0.5;

  float *A, *B, *C;
  float *A_blocked, *B_blocked, *C_blocked;

  struct multiply_stream_t *multiply_stream;

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
#endif

  /* Fill matrices with random data. */
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      A[get_index(N, i, j)] = rand()/(float) RAND_MAX;
      B[get_index(N, i, j)] = rand()/(float) RAND_MAX;
      C[get_index(N, i, j)] = rand()/(float) RAND_MAX;

      //A[get_index(N, i, j)] = get_index(N, i, j)+1;
      //B[get_index(N, i, j)] = get_index(N, i, j)+1;
      //C[get_index(N, i, j)] = get_index(N, i, j)+1;

      //printf("get_index(%u, %u, %u) = %u\n", N, i, j, get_index(N, i, j));
    }
  }

  //for (i = 0; i < N*N; i++)
  //{
  //  printf(" %i", (int) A[i]);
  //}
  //printf("\n");

  printf("A =\n");
  print_matrix(N, A);

  printf("B =\n");
  print_matrix(N, B);

  printf("C =\n");
  print_matrix(N, C);

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
  }

  stream_multiply(number_stream_elements, alpha, multiply_stream);

  printf("C =\n");
  print_matrix(N, C);

  /* Free memory. */
  free(A);
  free(B);
  free(C);
}
