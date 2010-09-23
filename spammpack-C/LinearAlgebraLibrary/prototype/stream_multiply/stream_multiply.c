#include "config.h"

#include <errno.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

#ifdef HAVE_XMMINTRIN_H
#include <xmmintrin.h>
#endif

#ifdef HAVE_PAPI
#include <papi.h>
#endif

#define N_BLOCK 4

//#undef HAVE_POSIX_MEMALIGN

//#define EXTERNAL_BLAS
//#define STREAM_KERNEL_1
//#define STREAM_KERNEL_2
//#define STREAM_KERNEL_3
//#define STREAM_KERNEL_4
//#define STREAM_KERNEL_5
//#define STREAM_KERNEL_6
//#define STREAM_KERNEL_7
//#define STREAM_KERNEL_8
//#define STREAM_KERNEL_9
//#define STREAM_KERNEL_10
//#define STREAM_KERNEL_11
//#define STREAM_KERNEL_12
//#define STREAM_KERNEL_13
//#define STREAM_KERNEL_14
//#define STREAM_KERNEL_15
//#define STREAM_KERNEL_16
//#define STREAM_KERNEL_17
//#define STREAM_KERNEL_18
//#define STREAM_KERNEL_19
#define STREAM_KERNEL_20
//#define POINTER_CHASE
//#define C_KERNEL
//#define NAIVE_KERNEL

#define STORE_DILATED_BLOCK

//#define DENSE_MULTIPLY

//#define DEBUG_LOOP

/* Allocate the dense matrix blocks in the tree in
   one contiguous chunk. The tree will then point into that chunk. */
#define CONTIGUOUS_ALLOCATE

struct matrix_box_t
{
  /* These boundaries are semi-open, i.e. the matrix indices covered in this
   * matrix node go from [M_lower, M_upper[ and [N_lower, N_upper[.
   */

  /* Row box. */
  unsigned int M_lower;
  unsigned int M_upper;

  /* Column box. */
  unsigned int N_lower;
  unsigned int N_upper;
};

struct matrix_node_t
{
  struct matrix_box_t box;
  struct matrix_node_t **child;
  float *block_dense;
#ifdef STORE_DILATED_BLOCK
  float *block_dilated;
#endif
};

struct matrix_t
{
  unsigned int N;
  struct matrix_box_t box;
  unsigned long long alignment;
  struct matrix_node_t *root;

#ifdef CONTIGUOUS_ALLOCATE
  float *contiguous;
#ifdef STORE_DILATED_BLOCK
  float *contiguous_dilated;
#endif
#endif
};

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
#if defined(STREAM_KERNEL_13) || defined(STREAM_KERNEL_14) || defined(STREAM_KERNEL_15) || defined(STREAM_KERNEL_16)
  char mask[8];
#elif defined(STREAM_KERNEL_17)
  float norm[8];
#elif defined(STREAM_KERNEL_18) || defined(STREAM_KERNEL_19) || defined(STREAM_KERNEL_20)
  float norm[32];
#endif
};

struct multiply_stream_index_t
{
  unsigned long long index;

  unsigned int A_index;
  unsigned int B_index;
  unsigned int C_index;
};

#ifdef EXTERNAL_BLAS
float
sgemm_ (char *, char *, int *, int *, int *, float *, float  *, int *,
    float *, int *, float *, float *, int *);
#endif

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

#ifdef STREAM_KERNEL_4
void
stream_kernel_4 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_5
void
stream_kernel_5 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_6
void
stream_kernel_6 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_7
void
stream_kernel_7 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_8
void
stream_kernel_8 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_9
void
stream_kernel_9 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_10
void
stream_kernel_10 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_11
void
stream_kernel_11 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_12
void
stream_kernel_12 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_13
void
stream_kernel_13 (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_14
void
stream_kernel_14_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_15
void
stream_kernel_15_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_16
void
stream_kernel_16_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_17
void
stream_kernel_17_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_18
void
stream_kernel_18_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_19
void
stream_kernel_19_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef STREAM_KERNEL_20
void
stream_kernel_20_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream);
#endif

#ifdef HAVE_PAPI
void
print_papi_error (const char *message, const int errorcode)
{
  char error[1001];

  PAPI_perror(errorcode, error, 1000);
  printf("PAPI error: %s - %s\n", message, error);
  exit(1);
}
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
    if (k & bitmask) { result |= addmask; }
    addmask <<= 1;
    if (j & bitmask) { result |= addmask; }
    addmask <<= 1;
    if (i & bitmask) { result |= addmask; }
    addmask <<= 1;
    bitmask <<= 1;
  }

  return result;
}

unsigned int
interleave_2_index (const unsigned int i, const unsigned int j)
{
  unsigned int result = 0;
  unsigned int bitmask = 1;
  unsigned int addmask = 1;
  unsigned int index;

  for (index = 0; index < sizeof(unsigned int)*8; index++)
  {
    if (j & bitmask) { result |= addmask; }
    addmask <<= 1;
    if (i & bitmask) { result |= addmask; }
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

unsigned int
get_Z_curve_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return interleave_2_index(i, j);
}

/* Returns the linear index of a submatrix block in a (N*N_BLOCK)x(N*N_BLOCK)
 * matrix. The blocks can be arranged in any order, the stream_kernel() does
 * not need to know what order this is.
 */
unsigned int
get_block_index (const unsigned int N, const unsigned int i, const unsigned int j) __attribute__ ((pure));

unsigned int
get_block_index (const unsigned int N, const unsigned int i, const unsigned int j)
{
  return get_Morton_index(N, i, j);
}

/* Returns the linear index of the submatrix block in a NxN matrix. */
unsigned int
get_blocked_index (const unsigned int N, const unsigned int i, const unsigned int j) __attribute__ ((pure));

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
get_index (const unsigned int N, const unsigned int i, const unsigned int j) __attribute__ ((pure));

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
    struct multiply_stream_t *multiply_stream,
    float **A_stream,
    float **B_stream,
    float **C_stream)
{
#if defined(STREAM_KERNEL_1)
  stream_kernel_1(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_2)
  stream_kernel_2(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_3)
  stream_kernel_3(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_4)
  stream_kernel_4(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_5)
  stream_kernel_5(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_6)
  stream_kernel_6(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_7)
  stream_kernel_7(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_8)
  stream_kernel_8(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_9)
  stream_kernel_9(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_10)
  stream_kernel_10(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_11)
  stream_kernel_11(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_12)
  stream_kernel_12(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_13)
  stream_kernel_13(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_14)
  stream_kernel_14_SSE_intrinsics(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_15)
  stream_kernel_15_SSE_intrinsics(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_16)
  stream_kernel_16_SSE_intrinsics(number_stream_elements, alpha, multiply_stream);

#elif defined(STREAM_KERNEL_17)
  stream_kernel_17_SSE_intrinsics(number_stream_elements, alpha, 1e-18, multiply_stream);

#elif defined(STREAM_KERNEL_18)
  stream_kernel_18_SSE_intrinsics(number_stream_elements, alpha, 1e-18, multiply_stream);

#elif defined(STREAM_KERNEL_19)
  stream_kernel_19_SSE_intrinsics(number_stream_elements, alpha, 1e-18, multiply_stream);

#elif defined(STREAM_KERNEL_20)
  stream_kernel_20_SSE_intrinsics(number_stream_elements, alpha, 1e-18, multiply_stream);

#elif defined(POINTER_CHASE)

#define READAHEAD 20

  unsigned long long stream_index;
  float *restrict A, *restrict B, *restrict C;

  if (helper_pid == 0)
  {
    /* The prefetch helper. */

  }

  else
  {
    /* Loops over all elements in stream. */
    for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
    {
      //if (stream_index < number_stream_elements-READAHEAD)
      //{
      //  _mm_prefetch(multiply_stream[stream_index+READAHEAD].A_block, _MM_HINT_T0);
      //  _mm_prefetch(multiply_stream[stream_index+READAHEAD].B_block, _MM_HINT_T0);
      //  _mm_prefetch(multiply_stream[stream_index+READAHEAD].C_block, _MM_HINT_T0);

      //  /* Get the helper to prefetch. */
      //}

      A = multiply_stream[stream_index].A_block;
      B = multiply_stream[stream_index].B_block;
      C = multiply_stream[stream_index].C_block;

      /* Do something so that the compiler does not optimize out the poitner
       * assignment and we are forcing the processor to load the matrix block.
       */
      if (A[0] == 0.0) { printf("do something\n"); }
      if (B[0] == 0.0) { printf("do something\n"); }
      if (C[0] == 0.0) { printf("do something\n"); }
    }
  }

#elif defined(C_KERNEL)

//#define READAHEAD 10
//#define C_KERNEL_VERSION_1
//#define C_KERNEL_VERSION_2
//#define C_KERNEL_VERSION_3
#define C_KERNEL_VERSION_4

#if defined(C_KERNEL_VERSION_1)
  short i;
  unsigned long long stream_index;
  float *restrict A, *restrict B, *restrict C;
  __m128 A_element, B_row, C_row, alpha_row;

  /* Load alpha. */
  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
  {
#ifdef READAHEAD
    if (stream_index < number_stream_elements-READAHEAD)
    {
      _mm_prefetch((void*) multiply_stream[stream_index+READAHEAD].A_block, _MM_HINT_T0);
      _mm_prefetch((void*) multiply_stream[stream_index+READAHEAD].B_block, _MM_HINT_T0);
      _mm_prefetch((void*) multiply_stream[stream_index+READAHEAD].C_block, _MM_HINT_T0);
    }
#endif

    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;

    for (i = 0; i < 4; i++)
    {
      A_element = _mm_load_ps(&A[(i*4+0)*4]);
      B_row = _mm_load_ps(&B[0*4]);
      C_row = _mm_mul_ps(A_element, B_row);

      A_element = _mm_load_ps(&A[(i*4+1)*4]);
      B_row = _mm_load_ps(&B[1*4]);
      C_row = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row);

      A_element = _mm_load_ps(&A[(i*4+2)*4]);
      B_row = _mm_load_ps(&B[2*4]);
      C_row = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row);

      A_element = _mm_load_ps(&A[(i*4+3)*4]);
      B_row = _mm_load_ps(&B[3*4]);
      C_row = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row);

      C_row = _mm_mul_ps(alpha_row, C_row);
      C_row = _mm_add_ps(_mm_load_ps(&C[i*4]), C_row);
      _mm_store_ps(&C[i*4], C_row);
    }
  }
#elif defined(C_KERNEL_VERSION_2)
  short i;
  unsigned long long stream_index;
  __m128 alpha_row;

  float *restrict A_1, *restrict B_1, *restrict C_1;
  __m128 A_element_1, B_row_1, C_row_1;

  float *restrict A_2, *restrict B_2, *restrict C_2;
  __m128 A_element_2, B_row_2, C_row_2;

  /* Load alpha. */
  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < number_stream_elements; stream_index += 2)
  {
    A_1 = multiply_stream[stream_index].A_block;
    B_1 = multiply_stream[stream_index].B_block;
    C_1 = multiply_stream[stream_index].C_block;

    for (i = 0; i < 4; i++)
    {
      A_element_1 = _mm_load_ps(&A_1[(i*4+0)*4]);
      B_row_1 = _mm_load_ps(&B_1[0*4]);
      C_row_1 = _mm_mul_ps(A_element_1, B_row_1);

      A_element_1 = _mm_load_ps(&A_1[(i*4+1)*4]);
      B_row_1 = _mm_load_ps(&B_1[1*4]);
      C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_1, B_row_1), C_row_1);

      A_element_1 = _mm_load_ps(&A_1[(i*4+2)*4]);
      B_row_1 = _mm_load_ps(&B_1[2*4]);
      C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_1, B_row_1), C_row_1);

      A_element_1 = _mm_load_ps(&A_1[(i*4+3)*4]);
      B_row_1 = _mm_load_ps(&B_1[3*4]);
      C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_1, B_row_1), C_row_1);

      C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
      C_row_1 = _mm_add_ps(_mm_load_ps(&C_1[i*4]), C_row_1);
      _mm_store_ps(&C_1[i*4], C_row_1);
    }

    A_2 = multiply_stream[stream_index+1].A_block;
    B_2 = multiply_stream[stream_index+1].B_block;
    C_2 = multiply_stream[stream_index+1].C_block;

    for (i = 0; i < 4; i++)
    {
      A_element_2 = _mm_load_ps(&A_2[(i*4+0)*4]);
      B_row_2 = _mm_load_ps(&B_2[0*4]);
      C_row_2 = _mm_mul_ps(A_element_2, B_row_2);

      A_element_2 = _mm_load_ps(&A_2[(i*4+1)*4]);
      B_row_2 = _mm_load_ps(&B_2[1*4]);
      C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_2, B_row_2), C_row_2);

      A_element_2 = _mm_load_ps(&A_2[(i*4+2)*4]);
      B_row_2 = _mm_load_ps(&B_2[2*4]);
      C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_2, B_row_2), C_row_2);

      A_element_2 = _mm_load_ps(&A_2[(i*4+3)*4]);
      B_row_2 = _mm_load_ps(&B_2[3*4]);
      C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_2, B_row_2), C_row_2);

      C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
      C_row_2 = _mm_add_ps(_mm_load_ps(&C_2[i*4]), C_row_2);
      _mm_store_ps(&C_2[i*4], C_row_2);
    }
  }
#elif defined(C_KERNEL_VERSION_3)
  float alpha_row[4] __attribute__((aligned(64)));

  /* Create 4-element float vector for alpha. */
  alpha_row[0] = alpha;
  alpha_row[1] = alpha;
  alpha_row[2] = alpha;
  alpha_row[3] = alpha;

  /* Store old rounding mode. */
  unsigned int old_xcsr = _mm_getcsr();

  /* Set DAZ and FTZ. */
  _mm_setcsr(old_xcsr | 0x8040);

  __asm__ volatile (

      /* Test whether to loop at all. */
      "cmp $0, %0\n\t"
      "jz done\n\t"

      /* Zero loop counter. */
      "xor %%rax, %%rax\n\t"

      /* Loop. */
      ".align 16\n"
      "loop:\n\t"

      /* Multiply stream block. */
      "mov 0x0(%2), %%r8\n\t"
      "mov 0x8(%2), %%r9\n\t"
      "mov 0x10(%2), %%r10\n\t"

      /* Load C matrix block. */
      "movaps (%%r10), %%xmm12\n\t"
      "movaps 4*4(%%r10), %%xmm13\n\t"
      "movaps 8*4(%%r10), %%xmm14\n\t"
      "movaps 16*4(%%r10), %%xmm15\n\t"

      /* Load A and B. */
      "movaps (%%r9), %%xmm0\n\t" /* Load B(1,1) B(1,2) B(1,3) B(1,4). */
      "movaps 4*4(%%r9), %%xmm1\n\t" /* Load B(2,1) B(2,2) B(2,3) B(2,4). */
      "movaps 8*4(%%r9), %%xmm2\n\t" /* Load B(3,1) B(3,2) B(3,3) B(3,4). */
      "movaps 16*4(%%r9), %%xmm3\n\t" /* Load B(4,1) B(4,2) B(4,3) B(4,4). */

      "movaps  0*4*4(%%r8), %%xmm4\n\t" /* Load A(1,1). */
      "movaps  1*4*4(%%r8), %%xmm5\n\t" /* Load A(1,2). */

      "movaps  2*4*4(%%r8), %%xmm6\n\t" /* Load A(1,3). */
      "movaps  3*4*4(%%r8), %%xmm7\n\t" /* Load A(1,4). */

      "movaps  4*4*4(%%r8), %%xmm8\n\t" /* Load A(2,1). */
      "movaps  5*4*4(%%r8), %%xmm9\n\t" /* Load A(2,2). */

      "movaps  6*4*4(%%r8), %%xmm10\n\t" /* Load A(2,3). */
      "movaps  7*4*4(%%r8), %%xmm11\n\t" /* Load A(2,4). */

      "mulps %%xmm0, %%xmm4\n\t"
      "mulps %%xmm1, %%xmm5\n\t"

      "mulps %%xmm3, %%xmm6\n\t"
      "mulps %%xmm4, %%xmm7\n\t"

      "mulps %%xmm0, %%xmm8\n\t"
      "mulps %%xmm1, %%xmm9\n\t"

      "mulps %%xmm3, %%xmm10\n\t"
      "mulps %%xmm4, %%xmm11\n\t"

      "addps %%xmm4, %%xmm12\n\t"
      "addps %%xmm5, %%xmm12\n\t"

      "addps %%xmm6, %%xmm13\n\t"
      "addps %%xmm7, %%xmm13\n\t"

      "addps %%xmm8, %%xmm14\n\t"
      "addps %%xmm9, %%xmm14\n\t"

      "addps %%xmm10, %%xmm15\n\t"
      "addps %%xmm11, %%xmm15\n\t"

      "movaps  8*4*4(%%r8), %%xmm4\n\t" /* Load A(3,1). */
      "movaps  9*4*4(%%r8), %%xmm5\n\t" /* Load A(3,2). */

      "movaps 10*4*4(%%r8), %%xmm6\n\t" /* Load A(3,3). */
      "movaps 11*4*4(%%r8), %%xmm7\n\t" /* Load A(3,4). */

      "movaps 12*4*4(%%r8), %%xmm8\n\t" /* Load A(4,1). */
      "movaps 13*4*4(%%r8), %%xmm9\n\t" /* Load A(4,2). */

      "movaps 14*4*4(%%r8), %%xmm10\n\t" /* Load A(4,3). */
      "movaps 15*4*4(%%r8), %%xmm11\n\t" /* Load A(4,4). */

      "mulps %%xmm0, %%xmm4\n\t"
      "mulps %%xmm1, %%xmm5\n\t"

      "mulps %%xmm3, %%xmm6\n\t"
      "mulps %%xmm4, %%xmm7\n\t"

      "mulps %%xmm0, %%xmm8\n\t"
      "mulps %%xmm1, %%xmm9\n\t"

      "mulps %%xmm3, %%xmm10\n\t"
      "mulps %%xmm4, %%xmm11\n\t"

      "addps %%xmm4, %%xmm12\n\t"
      "addps %%xmm5, %%xmm12\n\t"

      "addps %%xmm6, %%xmm13\n\t"
      "addps %%xmm7, %%xmm13\n\t"

      "addps %%xmm8, %%xmm14\n\t"
      "addps %%xmm9, %%xmm14\n\t"

      "addps %%xmm10, %%xmm15\n\t"
      "addps %%xmm11, %%xmm15\n\t"

      /* Write out C. */
      "movaps %%xmm12, (%%r10)\n\t"
      "movaps %%xmm13, 4*4(%%r10)\n\t"
      "movaps %%xmm14, 8*4(%%r10)\n\t"
      "movaps %%xmm15, 16*4(%%r10)\n\t"

      /* Increment loop counter. */
      "add $1, %%rax\n\t"
      "cmp %0, %%rax\n\t"
      "jl loop\n\t"

      /* Done. */
      ".align 16\n"
      "done:\n"

      : /* output registers. */
      : "r" (number_stream_elements),
        "r" (alpha_row),
        "r" (multiply_stream) /* Input registers. */
      : "memory", "rax", "r8", "r9", "r10"
      );

  /* Restore old rounding mode. */
  _mm_setcsr(old_xcsr);
#elif defined(C_KERNEL_VERSION_4)
  short i;
  unsigned long long stream_index;
  float *restrict A, *restrict B, *restrict C;

  __m128 A_element_11, A_element_12, A_element_13, A_element_14;
  __m128 A_element_21, A_element_22, A_element_23, A_element_24;
  __m128 A_element_31, A_element_32, A_element_33, A_element_34;
  __m128 A_element_41, A_element_42, A_element_43, A_element_44;

  __m128 B_row_1, B_row_2, B_row_3, B_row_4;
  __m128 C_row_1, C_row_2, C_row_3, C_row_4;

  __m128 alpha_row;

  /* Load alpha. */
  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < number_stream_elements; stream_index++)
  {
#ifdef READAHEAD
    if (stream_index < number_stream_elements-READAHEAD)
    {
      _mm_prefetch((void*) multiply_stream[stream_index+READAHEAD].A_block, _MM_HINT_T0);
      _mm_prefetch((void*) multiply_stream[stream_index+READAHEAD].B_block, _MM_HINT_T0);
      _mm_prefetch((void*) multiply_stream[stream_index+READAHEAD].C_block, _MM_HINT_T0);
    }
#endif

    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;

    A_element_11 = _mm_load_ps(&A[(0*4+0)*4]);
    A_element_12 = _mm_load_ps(&A[(0*4+1)*4]);
    A_element_13 = _mm_load_ps(&A[(0*4+2)*4]);
    A_element_14 = _mm_load_ps(&A[(0*4+3)*4]);

    A_element_21 = _mm_load_ps(&A[(1*4+0)*4]);
    A_element_22 = _mm_load_ps(&A[(1*4+1)*4]);
    A_element_23 = _mm_load_ps(&A[(1*4+2)*4]);
    A_element_24 = _mm_load_ps(&A[(1*4+3)*4]);

    A_element_31 = _mm_load_ps(&A[(2*4+0)*4]);
    A_element_32 = _mm_load_ps(&A[(2*4+1)*4]);
    A_element_33 = _mm_load_ps(&A[(2*4+2)*4]);
    A_element_34 = _mm_load_ps(&A[(2*4+3)*4]);

    A_element_41 = _mm_load_ps(&A[(3*4+0)*4]);
    A_element_42 = _mm_load_ps(&A[(3*4+1)*4]);
    A_element_43 = _mm_load_ps(&A[(3*4+2)*4]);
    A_element_44 = _mm_load_ps(&A[(3*4+3)*4]);

    B_row_1 = _mm_load_ps(&B[0*4]);
    B_row_2 = _mm_load_ps(&B[1*4]);
    B_row_3 = _mm_load_ps(&B[2*4]);
    B_row_4 = _mm_load_ps(&B[3*4]);

    C_row_1 = _mm_mul_ps(A_element_11, B_row_1);
    C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_12, B_row_2), C_row_1);
    C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_13, B_row_3), C_row_1);
    C_row_1 = _mm_add_ps(_mm_mul_ps(A_element_14, B_row_4), C_row_1);
    C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
    C_row_1 = _mm_add_ps(_mm_load_ps(&C[0*4]), C_row_1);
    _mm_store_ps(&C[0*4], C_row_1);

    C_row_2 = _mm_mul_ps(A_element_21, B_row_1);
    C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_22, B_row_2), C_row_2);
    C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_23, B_row_3), C_row_2);
    C_row_2 = _mm_add_ps(_mm_mul_ps(A_element_24, B_row_4), C_row_2);
    C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
    C_row_2 = _mm_add_ps(_mm_load_ps(&C[1*4]), C_row_2);
    _mm_store_ps(&C[1*4], C_row_2);

    C_row_3 = _mm_mul_ps(A_element_31, B_row_1);
    C_row_3 = _mm_add_ps(_mm_mul_ps(A_element_32, B_row_2), C_row_3);
    C_row_3 = _mm_add_ps(_mm_mul_ps(A_element_33, B_row_3), C_row_3);
    C_row_3 = _mm_add_ps(_mm_mul_ps(A_element_34, B_row_4), C_row_3);
    C_row_3 = _mm_mul_ps(alpha_row, C_row_3);
    C_row_3 = _mm_add_ps(_mm_load_ps(&C[2*4]), C_row_3);
    _mm_store_ps(&C[2*4], C_row_3);

    C_row_4 = _mm_mul_ps(A_element_41, B_row_1);
    C_row_4 = _mm_add_ps(_mm_mul_ps(A_element_42, B_row_2), C_row_4);
    C_row_4 = _mm_add_ps(_mm_mul_ps(A_element_43, B_row_3), C_row_4);
    C_row_4 = _mm_add_ps(_mm_mul_ps(A_element_44, B_row_4), C_row_4);
    C_row_4 = _mm_mul_ps(alpha_row, C_row_4);
    C_row_4 = _mm_add_ps(_mm_load_ps(&C[3*4]), C_row_4);
    _mm_store_ps(&C[3*4], C_row_4);
  }
#endif

#elif defined(NAIVE_KERNEL)
  unsigned int i, j, k;

  /* Multiply stream. */
  for (stream_index = 0; stream_index < number_stream_elements; stream_index++) {
    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;

    for (i = 0; i < N_BLOCK; i++) {
      for (j = 0; j < N_BLOCK; j++) {
        for (k = 0; k < N_BLOCK; k++)
        {
          C[get_block_index(N_BLOCK, i, j)] += alpha*A[get_block_index(N_BLOCK, i, k)]*B[get_block_index(N_BLOCK, k, j)];
        }
      }
    }
  }

#else
  printf("no kernel\n");
  exit(1);

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

struct matrix_t*
spamm_new (const unsigned int N, const unsigned int alignment) __attribute__ ((malloc));

struct matrix_t*
spamm_new (const unsigned int N, const unsigned int alignment)
{
#ifdef HAVE_POSIX_MEMALIGN
  int allocation_result;
#endif
  struct matrix_t *new_matrix;

  if ((new_matrix = (struct matrix_t*) malloc(sizeof(struct matrix_t))) == NULL)
  {
    printf("error allocating new matrix\n");
    exit(1);
  }

  new_matrix->N = N;
  new_matrix->box.M_lower = 0;
  new_matrix->box.N_lower = 0;
  new_matrix->box.M_upper = N;
  new_matrix->box.N_upper = N;
  new_matrix->alignment = alignment;
  new_matrix->root = NULL;

#ifdef CONTIGUOUS_ALLOCATE

#ifdef HAVE_POSIX_MEMALIGN
  if ((allocation_result = posix_memalign((void**) &new_matrix->contiguous, alignment, sizeof(float)*N*N)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

#else
  if ((new_matrix->contiguous = (float*) malloc(sizeof(float)*N*N)) == NULL)
  {
    printf("error allocating contiguous memory for matrix\n");
    exit(1);
  }

#endif

#ifdef STORE_DILATED_BLOCK

#ifdef HAVE_POSIX_MEMALIGN
  if ((allocation_result = posix_memalign((void**) &new_matrix->contiguous_dilated, alignment, sizeof(float)*N*N*4)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

#else
  if ((new_matrix->contiguous_dilated = (float*) malloc(sizeof(float)*N*N*4)) == NULL)
  {
    printf("error allocating contiguous_dilated memory for matrix\n");
    exit(1);
  }

#endif

#endif

#endif

  return new_matrix;
}

struct matrix_node_t*
spamm_new_node (const unsigned int M_lower, const unsigned int M_upper,
    const unsigned int N_lower, const unsigned int N_upper,
    const unsigned long long alignment) __attribute__ ((malloc));

struct matrix_node_t*
spamm_new_node (const unsigned int M_lower, const unsigned int M_upper,
    const unsigned int N_lower, const unsigned int N_upper,
    const unsigned long long alignment)
{
  struct matrix_node_t *new_node;

  if ((new_node = (struct matrix_node_t*) malloc(sizeof(struct matrix_node_t))) == NULL)
  {
    printf("error allocating new node\n");
    exit(1);
  }

  new_node->box.M_lower = M_lower;
  new_node->box.M_upper = M_upper;
  new_node->box.N_lower = N_lower;
  new_node->box.N_upper = N_upper;
  new_node->child = NULL;
  new_node->block_dense = NULL;
#ifdef STORE_DILATED_BLOCK
  new_node->block_dilated = NULL;
#endif

  return new_node;
}

void
spamm_set_node (struct matrix_node_t *node, struct matrix_t *A,
    const unsigned int i, const unsigned int j, const float Aij,
    const unsigned long long alignment)
{
  short int row_index, column_index;
  short int child_index;

  if (node->box.M_upper-node->box.M_lower == N_BLOCK &&
      node->box.N_upper-node->box.N_lower == N_BLOCK)
  {
#ifdef CONTIGUOUS_ALLOCATE

    if (node->block_dense == NULL)
    {
      /* Figure out where to point into the contiguous region. */

      node->block_dense = &A->contiguous[interleave_2_index(i/N_BLOCK, j/N_BLOCK)*N_BLOCK*N_BLOCK];
#ifdef STORE_DILATED_BLOCK
      node->block_dilated = &A->contiguous_dilated[interleave_2_index(i/N_BLOCK, j/N_BLOCK)*N_BLOCK*N_BLOCK*4];
#endif

#else
    if (node->block_dense == NULL)
    {
      /* Allocate new dense matrix block. */
#if defined(HAVE_POSIX_MEMALIGN)
      int allocation_result;
      if ((allocation_result = posix_memalign((void**) &node->block_dense, alignment, sizeof(float)*N_BLOCK*N_BLOCK)) != 0)
      {
        switch (allocation_result)
        {
          case EINVAL:
            printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
            break;

          case ENOMEM:
            printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
            break;

          default:
            printf("Unknown error code\n");
            break;
        }
        exit(1);
      }

#ifdef STORE_DILATED_BLOCK
      if ((allocation_result = posix_memalign((void**) &node->block_dilated, alignment, sizeof(float)*4*N_BLOCK*N_BLOCK)) != 0)
      {
        switch (allocation_result)
        {
          case EINVAL:
            printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
            break;

          case ENOMEM:
            printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
            break;

          default:
            printf("Unknown error code\n");
            break;
        }
        exit(1);
      }
#endif

#else
      if ((node->block_dense = (float*) malloc(sizeof(float)*N_BLOCK*N_BLOCK)) == NULL)
      {
        printf("error allocating new dense matrix block\n");
        exit(1);
      }

#ifdef STORE_DILATED_BLOCK
      if ((node->block_dilated = (float*) malloc(sizeof(float)*4*N_BLOCK*N_BLOCK)) == NULL)
      {
        printf("error allocating new dense matrix block\n");
        exit(1);
      }
#endif

#endif

#endif

      /* Zero out the newly allocated block. */
      for (row_index = 0; row_index < N_BLOCK; row_index++) {
        for (column_index = 0; column_index < N_BLOCK; column_index++)
        {
          node->block_dense[row_index*N_BLOCK+column_index] = 0.0;
#ifdef STORE_DILATED_BLOCK
          node->block_dilated[(row_index*N_BLOCK+column_index)*4+0] = 0.0;
          node->block_dilated[(row_index*N_BLOCK+column_index)*4+1] = 0.0;
          node->block_dilated[(row_index*N_BLOCK+column_index)*4+2] = 0.0;
          node->block_dilated[(row_index*N_BLOCK+column_index)*4+3] = 0.0;
#endif
        }
      }
    }

    node->block_dense[(i-node->box.M_lower)*N_BLOCK+(j-node->box.N_lower)] = Aij;
#ifdef STORE_DILATED_BLOCK
    node->block_dilated[((i-node->box.M_lower)*N_BLOCK+(j-node->box.N_lower))*4+0] = Aij;
    node->block_dilated[((i-node->box.M_lower)*N_BLOCK+(j-node->box.N_lower))*4+1] = Aij;
    node->block_dilated[((i-node->box.M_lower)*N_BLOCK+(j-node->box.N_lower))*4+2] = Aij;
    node->block_dilated[((i-node->box.M_lower)*N_BLOCK+(j-node->box.N_lower))*4+3] = Aij;
#endif
  }

  else
  {
    if (node->child == NULL)
    {
      if ((node->child = (struct matrix_node_t**) malloc(sizeof(struct matrix_node_t*)*4)) == NULL)
      {
        printf("error allocating children pointers\n");
        exit(1);
      }

      node->child[0*2+0] = spamm_new_node(node->box.M_lower,
          node->box.M_lower+(node->box.M_upper-node->box.M_lower)/2,
          node->box.N_lower,
          node->box.N_lower+(node->box.N_upper-node->box.N_lower)/2,
          alignment);
      node->child[0*2+1] = spamm_new_node(node->box.M_lower,
          node->box.M_lower+(node->box.M_upper-node->box.M_lower)/2,
          node->box.N_lower+(node->box.N_upper-node->box.N_lower)/2,
          node->box.N_upper,
          alignment);
      node->child[1*2+0] = spamm_new_node(node->box.M_lower+(node->box.M_upper-node->box.M_lower)/2,
          node->box.M_upper,
          node->box.N_lower,
          node->box.N_lower+(node->box.N_upper-node->box.N_lower)/2,
          alignment);
      node->child[1*2+1] = spamm_new_node(node->box.M_lower+(node->box.M_upper-node->box.M_lower)/2,
          node->box.M_upper,
          node->box.N_lower+(node->box.N_upper-node->box.N_lower)/2,
          node->box.N_upper,
          alignment);
    }

    if (i < node->box.M_lower+(node->box.M_upper-node->box.M_lower)/2)
    {
      child_index = 0;
    }
    else
    {
      child_index = 2;
    }

    if (j >= node->box.N_lower+(node->box.N_upper-node->box.N_lower)/2)
    {
      child_index |= 1;
    }

    /* Descent to child node. */
    spamm_set_node(node->child[child_index], A, i, j, Aij, alignment);
  }
}

void
spamm_set (struct matrix_t *A, const unsigned int i, const unsigned int j, const float Aij)
{
  if (A->root == NULL)
  {
    A->root = spamm_new_node(A->box.M_lower, A->box.M_upper, A->box.N_lower, A->box.N_upper, A->alignment);
  }

  spamm_set_node(A->root, A, i, j, Aij, A->alignment);
}

float
spamm_get_node (struct matrix_node_t *node, const unsigned int i, const unsigned int j)
{
  short int child_index;

  if (node->child != NULL)
  {
    if (i < node->box.M_lower+(node->box.M_upper-node->box.M_lower)/2)
    {
      child_index = 0;
    }
    else
    {
      child_index = 2;
    }

    if (j >= node->box.N_lower+(node->box.N_upper-node->box.N_lower)/2)
    {
      child_index |= 1;
    }

    /* Descent to child node. */
    return spamm_get_node(node->child[child_index], i, j);
  }

  else if (node->block_dense != NULL)
  {
    return node->block_dense[(i-node->box.M_lower)*N_BLOCK+(j-node->box.N_lower)];
  }

  else { return 0.0; }
}

float
spamm_get (const struct matrix_t *A, const unsigned int i, const unsigned int j)
{
  if (A->root == NULL) { return 0; }
  else
  {
    return spamm_get_node(A->root, i, j);
  }
}

void
spamm_print (struct matrix_t *A)
{
  unsigned int i, j;

  for (i = 0; i < A->N; i++) {
    for (j = 0; j < A->N; j++)
    {
      printf(" %f", spamm_get(A, i, j));
    }
    printf("\n");
  }

#ifdef CONTIGUOUS_ALLOCATE
  printf("linear: ");
  for (i = 0; i < A->N*A->N; i++)
  {
    printf(" %f", A->contiguous[i]);
  }
  printf("\n");
#endif
}

void
spamm_multiply_node (const struct matrix_node_t *A_node,
    const struct matrix_node_t *B_node,
    const struct matrix_node_t *C_node,
    unsigned long long *index,
    struct multiply_stream_t *stream,
    struct multiply_stream_index_t *stream_index,
    float **A_stream,
    float **B_stream,
    float **C_stream)
{
  int i, j, k;

  if (A_node->block_dense != NULL)
  {
#ifdef STORE_DILATED_BLOCK
    stream[*index].A_block = A_node->block_dilated;
#else
    stream[*index].A_block = A_node->block_dense;
#endif
    stream[*index].B_block = B_node->block_dense;
    stream[*index].C_block = C_node->block_dense;

#ifdef STORE_DILATED_BLOCK
    A_stream[*index] = A_node->block_dilated;
#else
    A_stream[*index] = A_node->block_dense;
#endif
    B_stream[*index] = B_node->block_dense;
    C_stream[*index] = C_node->block_dense;

#if defined (STREAM_KERNEL_13) || defined(STREAM_KERNEL_14) || defined(STREAM_KERNEL_15) || defined(STREAM_KERNEL_16)
    stream[*index].mask[0] = 1;
    stream[*index].mask[1] = 1;
    stream[*index].mask[2] = 1;
    stream[*index].mask[3] = 1;
    stream[*index].mask[4] = 1;
    stream[*index].mask[5] = 1;
    stream[*index].mask[6] = 1;
    stream[*index].mask[7] = 1;
#elif defined(STREAM_KERNEL_17)
    for (i = 0; i < 8; i++)
    {
      stream[*index].norm[i] = rand()/(float) RAND_MAX;
    }
#elif defined(STREAM_KERNEL_18) || defined(STREAM_KERNEL_19) || defined(STREAM_KERNEL_20)
    for (i = 0; i < 32; i++)
    {
      stream[*index].norm[i] = rand()/(float) RAND_MAX;
    }
#endif

    (*index)++;
  }

  else
  {
    for (i = 0; i < 2; i++) {
      for (j = 0; j < 2; j++) {
        for (k = 0; k < 2; k++)
        {
          spamm_multiply_node(A_node->child[i*2+k],
              B_node->child[k*2+j],
              C_node->child[i*2+j],
              index, stream, stream_index, A_stream, B_stream, C_stream);
        }
      }
    }
  }
}

void
spamm_multiply (const struct matrix_t *A,
    const struct matrix_t *B,
    const struct matrix_t *C,
    unsigned long long *index,
    struct multiply_stream_t *stream,
    struct multiply_stream_index_t *stream_index,
    float **A_stream,
    float **B_stream,
    float **C_stream)
{
  /* Recursively multiply. */
  spamm_multiply_node(A->root, B->root, C->root, index, stream, stream_index, A_stream, B_stream, C_stream);
}

float
spamm_compare_matrix (const struct matrix_t *A, float *B)
{
  unsigned int i, j;
  float max_diff = 0;

  for (i = 0; i < A->N; i++) {
    for (j = 0; j < A->N; j++)
    {
      if (fabs(spamm_get(A, i, j)-B[get_index(A->N, i, j)]) > max_diff)
      {
        max_diff = fabs(spamm_get(A, i, j)-B[get_index(A->N, i, j)]);
      }
    }
  }

  return max_diff;
}

void
spamm_free_node (struct matrix_node_t **node)
{
  short int i;

  if ((*node)->child != NULL)
  {
    for (i = 0; i < 4; i++)
    {
      spamm_free_node(&(*node)->child[i]);
    }
  }

  else if ((*node)->block_dense != NULL)
  {
#if ! defined(CONTIGUOUS_ALLOCATE)
    free((*node)->block_dense);
#endif
    (*node)->block_dense = NULL;

#ifdef STORE_DILATED_BLOCK
#if ! defined(CONTIGUOUS_ALLOCATE)
    free((*node)->block_dilated);
#endif
    (*node)->block_dilated = NULL;
#endif
  }

  free(*node);
  *node = NULL;
}

void
spamm_free (struct matrix_t **A)
{
  if ((*A)->root != NULL)
  {
    spamm_free_node(&(*A)->root);
  }
#ifdef CONTIGUOUS_ALLOCATE
  free((*A)->contiguous);
#ifdef STORE_DILATED_BLOCK
  free((*A)->contiguous_dilated);
#endif
#endif
  free(*A);
  *A = NULL;
}

void
print_usage (unsigned int N, unsigned int alignment, unsigned int long loops,
    unsigned int N_only_A, unsigned int N_only_B, unsigned int N_only_C)
{
  printf("Usage:\n");
  printf("\n");
  printf("-h            This help\n");
  printf("-N N          Use NxN matrices (default N = %u)\n", N);
  printf("--align N     Align memory buffer on N byte boundary (default N = %u)\n", alignment);
  printf("--loops N     Repeat each access test N times (default N = %llu)\n", loops);
  printf("--verify      Verify result\n");
  printf("--print       Print matrices\n");
  printf("--no-random   Full matrices with index values as opposed to random\n");
  printf("--sort        Sort stream\n");
  printf("--only_A N    Map only N matrix blocks in stream A (default N = %u)\n", N_only_A);
  printf("--only_B N    Map only N matrix blocks in stream A (default N = %u)\n", N_only_B);
  printf("--only_C N    Map only N matrix blocks in stream A (default N = %u)\n", N_only_C);
  printf("--histogram   Calculate address distance histogram of blocks in stream\n");
  printf("--seed N      Use N as the seed for rand() (default N = time())\n");
#ifdef HAVE_PAPI
  printf("--TOT_INS     Measure total instructions\n");
  printf("--TOT_CYC     Measure total cycles\n");
  printf("--RES_STL     Measure stalled cycles\n");
  printf("--L1_ICM      Measure L1 instruction misses\n");
  printf("--L1_DCM      Measure L1 data misses\n");
  printf("--L1_LDM      Measure L1 data load misses\n");
  printf("--L1_STM      Measure L1 data store misses\n");
  printf("--L1_DCH      Measure L1 data hits\n");
  printf("--L1_DCA      Measure L1 data accesses\n");
  printf("--L2_ICM      Measure L2 instruction misses\n");
  printf("--L2_DCM      Measure L2 data misses\n");
  printf("--L2_LDM      Measure L2 data load misses\n");
  printf("--L2_STM      Measure L2 data store misses\n");
  printf("--L2_DCH      Measure L2 data hits\n");
  printf("--L2_DCA      Measure L2 data accesses\n");
  printf("--TLB_IM      Measure TLB instruction misses\n");
  printf("--TLB_DM      Measure TLB data misses\n");
  printf("--TLB_SD      Measure TLB shootdowns\n");
#endif
}

int
main (int argc, char **argv)
{
  unsigned int N = 1024;
  unsigned int N_padded = N;

  unsigned int N_only_A = N*N;
  unsigned int N_only_B = N*N;
  unsigned int N_only_C = N*N;

  unsigned int alignment = 64;

#ifdef HAVE_POSIX_MEMALIGN
  int allocation_result;
#endif

  unsigned int i, j, k;
  unsigned int (*index_pairs)[2];
  unsigned long long stream_index;
  unsigned long long stream_index_A, stream_index_B, stream_index_C;
  unsigned long long number_stream_elements;

  float alpha = 1.2;
  float beta = 0.5;

  unsigned int rand_seed = 1;

  int verify = 0;
  double max_diff;

  int print = 0;

  int random_elements = 1;

  int sort_stream = 0;

  int histogram = 0;

  unsigned int loop;
  unsigned long long loops = 1;
  struct timeval start, stop;
  struct rusage rusage_start, rusage_stop;
  double walltime, usertime, systime, flops;

#ifdef HAVE_PAPI
  int papi_events = PAPI_NULL;
  long long *papi_values;
  int papi_value_index;
  int papi_errorcode;

  int load_TOT_INS = 0;
  int load_TOT_CYC = 0;
  int load_RES_STL = 0;
  int load_L1_ICM = 0;
  int load_L1_DCM = 0;
  int load_L1_LDM = 0;
  int load_L1_STM = 0;
  int load_L1_DCH = 0;
  int load_L1_DCA = 0;
  int load_L2_ICM = 0;
  int load_L2_DCM = 0;
  int load_L2_LDM = 0;
  int load_L2_STM = 0;
  int load_L2_DCH = 0;
  int load_L2_DCA = 0;
  int load_TLB_IM = 0;
  int load_TLB_DM = 0;
  int load_TLB_SD = 0;
#endif

  struct matrix_t *A_spamm, *B_spamm, *C_spamm;

  float *A, *B, *C, *D;

  float **A_stream, **B_stream, **C_stream;
  struct multiply_stream_t *multiply_stream;
  struct multiply_stream_index_t *multiply_stream_index;

  unsigned int unique_block_index;
  void **unique_blocks;
  unsigned long long *histogram_distances;
  unsigned int *histogram_counts;

  unsigned long long tree_size;
  unsigned int tree_depth;

  void *temp_pointer;
  unsigned int temp_unsigned_int;

  unsigned int memory_access_length, old_A_index, old_B_index, old_C_index;

  int parse;
  int longindex;
  char *short_options = "hN:";
  struct option long_options[] = {
    { "help", no_argument, NULL, 'h' },
    { "N", required_argument, NULL, 'N' },
    { "loops", required_argument, NULL, 'l' },
    { "align", required_argument, NULL, 'a' },
    { "verify", no_argument, NULL, 'v' },
    { "print", no_argument, NULL, 'p' },
    { "no-random", no_argument, NULL, 'r' },
    { "sort", no_argument, NULL, 's' },
    { "only_A", required_argument, NULL, 'x' },
    { "only_B", required_argument, NULL, 'y' },
    { "only_C", required_argument, NULL, 'z' },
    { "histogram", no_argument, NULL, 'q' },
    { "seed", required_argument, NULL, 'S' },
#ifdef HAVE_PAPI
    { "TOT_INS",  no_argument, NULL, '1' },
    { "TOT_CYC",  no_argument, NULL, '2' },
    { "RES_STL",  no_argument, NULL, '3' },
    { "L1_ICM",   no_argument, NULL, '4' },
    { "L1_DCM",   no_argument, NULL, '5' },
    { "L1_LDM",   no_argument, NULL, '6' },
    { "L1_STM",   no_argument, NULL, '7' },
    { "L1_DCH",   no_argument, NULL, '8' },
    { "L1_DCA",   no_argument, NULL, '9' },
    { "L2_ICM",   no_argument, NULL, 'A' },
    { "L2_DCM",   no_argument, NULL, 'B' },
    { "L2_LDM",   no_argument, NULL, 'C' },
    { "L2_STM",   no_argument, NULL, 'D' },
    { "L2_DCH",   no_argument, NULL, 'E' },
    { "L2_DCA",   no_argument, NULL, 'F' },
    { "TLB_IM",   no_argument, NULL, 'G' },
    { "TLB_DM",   no_argument, NULL, 'H' },
    { "TLB_SD",   no_argument, NULL, 'I' },
#endif
    { NULL, 0, NULL, 0 }
  };

  /* Read command line. */
  while ((parse = getopt_long(argc, argv, short_options, long_options, &longindex)) != -1)
  {
    switch (parse)
    {
      case 'h':
        print_usage(N, alignment, loops, N_only_A, N_only_B, N_only_C);
        return 0;
        break;

      case 'N':
        N = strtol(optarg, NULL, 10);
        break;

      case 'l':
        loops = strtol(optarg, NULL, 10);
        break;

      case 'a':
        alignment = strtol(optarg, NULL, 10);
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

      case 'x':
        N_only_A = strtol(optarg, NULL, 10);
        break;

      case 'y':
        N_only_B = strtol(optarg, NULL, 10);
        break;

      case 'z':
        N_only_C = strtol(optarg, NULL, 10);
        break;

      case 'q':
        histogram = 1;
        break;

      case 'S':
        rand_seed = strtol(optarg, NULL, 10);
        break;

#ifdef HAVE_PAPI
      case '1':
        load_TOT_INS = 1;
        break;

      case '2':
        load_TOT_CYC = 1;
        break;

      case '3':
        load_RES_STL = 1;
        break;

      case '4':
        load_L1_ICM = 1;
        break;

      case '5':
        load_L1_DCM = 1;
        break;

      case '6':
        load_L1_LDM = 1;
        break;

      case '7':
        load_L1_STM = 1;
        break;

      case '8':
        load_L1_DCH = 1;
        break;

      case '9':
        load_L1_DCA = 1;
        break;

      case 'A':
        load_L2_ICM = 1;
        break;

      case 'B':
        load_L2_DCM = 1;
        break;

      case 'C':
        load_L2_LDM = 1;
        break;

      case 'D':
        load_L2_STM = 1;
        break;

      case 'E':
        load_L2_DCH = 1;
        break;

      case 'F':
        load_L2_DCA = 1;
        break;

      case 'G':
        load_TLB_IM = 1;
        break;

      case 'H':
        load_TLB_DM = 1;
        break;

      case 'I':
        load_TLB_SD = 1;
        break;
#endif

      default:
        printf("unknown command line argument\n");
        print_usage(N, alignment, loops, N_only_A, N_only_B, N_only_C);
        return -1;
        break;
    }
  }

  if (optind < argc) {
    printf("non-option command line input: ");
    while (optind < argc)
    {
      printf("%s ", argv[optind++]);
    }
    printf("\n");
    exit(1);
  }

#ifdef HAVE_PAPI
  /* Do some PAPI. */
  if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
  {
    printf("can not initialize PAPI\n");
    exit(1);
  }

  papi_values = (long long*) malloc(sizeof(long long)*18);
  if ((papi_errorcode = PAPI_set_granularity(PAPI_GRN_MIN) != PAPI_OK)) { print_papi_error("failed to set granularity",  papi_errorcode); }
  if ((papi_errorcode = PAPI_create_eventset(&papi_events) != PAPI_OK)) { print_papi_error("failed to create eventset",  papi_errorcode); }
  if ((papi_errorcode = PAPI_assign_eventset_component(papi_events, 0) != PAPI_OK)) { print_papi_error("failed to assign component to eventset",  papi_errorcode); }

  if (load_TOT_INS && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TOT_INS))) { print_papi_error("failed to add PAPI_TOT_INS", papi_errorcode); }
  if (load_TOT_CYC && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TOT_CYC))) { print_papi_error("failed to add PAPI_TOT_CYC", papi_errorcode); }
  if (load_RES_STL && (papi_errorcode = PAPI_add_event(papi_events, PAPI_RES_STL))) { print_papi_error("failed to add PAPI_RES_STL", papi_errorcode); }

  if (load_L1_ICM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_ICM))) { print_papi_error("failed to add PAPI_L1_ICM",  papi_errorcode); }
  if (load_L1_DCM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_DCM))) { print_papi_error("failed to add PAPI_L1_DCM",  papi_errorcode); }
  if (load_L1_LDM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_LDM))) { print_papi_error("failed to add PAPI_L1_LDM",  papi_errorcode); }
  if (load_L1_STM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_STM))) { print_papi_error("failed to add PAPI_L1_STM",  papi_errorcode); }
  if (load_L1_DCH && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_DCH))) { print_papi_error("failed to add PAPI_L1_DCH",  papi_errorcode); }
  if (load_L1_DCA && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L1_DCA))) { print_papi_error("failed to add PAPI_L1_DCA",  papi_errorcode); }

  if (load_L2_ICM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_ICM))) { print_papi_error("failed to add PAPI_L2_ICM",  papi_errorcode); }
  if (load_L2_DCM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_DCM))) { print_papi_error("failed to add PAPI_L2_DCM",  papi_errorcode); }
  if (load_L2_LDM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_LDM))) { print_papi_error("failed to add PAPI_L2_LDM",  papi_errorcode); }
  if (load_L2_STM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_STM))) { print_papi_error("failed to add PAPI_L2_STM",  papi_errorcode); }
  if (load_L2_DCH && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_DCH))) { print_papi_error("failed to add PAPI_L2_DCH",  papi_errorcode); }
  if (load_L2_DCA && (papi_errorcode = PAPI_add_event(papi_events, PAPI_L2_DCA))) { print_papi_error("failed to add PAPI_L2_DCA",  papi_errorcode); }

  if (load_TLB_IM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TLB_IM))) { print_papi_error("failed to add PAPI_TLB_IM",  papi_errorcode); }
  if (load_TLB_DM && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TLB_DM))) { print_papi_error("failed to add PAPI_TLB_DM",  papi_errorcode); }
  if (load_TLB_SD && (papi_errorcode = PAPI_add_event(papi_events, PAPI_TLB_SD))) { print_papi_error("failed to add PAPI_TLB_SD",  papi_errorcode); }
#endif

#if defined(EXTERNAL_BLAS)
  printf("external blas\n");

#elif defined(STREAM_KERNEL_1)
  printf("using stream_kernel_1\n");

#elif defined(STREAM_KERNEL_2)
  printf("using stream_kernel_2\n");

#elif defined(STREAM_KERNEL_3)
  printf("using stream_kernel_3\n");

#elif defined(STREAM_KERNEL_4)
  printf("using stream_kernel_4\n");

#elif defined(STREAM_KERNEL_5)
  printf("using stream_kernel_5\n");

#elif defined(STREAM_KERNEL_6)
  printf("using stream_kernel_6\n");

#elif defined(STREAM_KERNEL_7)
  printf("using stream_kernel_7\n");

#elif defined(STREAM_KERNEL_8)
  printf("using stream_kernel_8\n");

#elif defined(STREAM_KERNEL_9)
  printf("using stream_kernel_9\n");

#elif defined(STREAM_KERNEL_10)
  printf("using stream_kernel_10\n");

#elif defined(STREAM_KERNEL_11)
  printf("using stream_kernel_11\n");

#elif defined(STREAM_KERNEL_12)
  printf("using stream_kernel_12\n");

#elif defined(STREAM_KERNEL_13)
  printf("using stream_kernel_13\n");

#elif defined(STREAM_KERNEL_14)
  printf("using stream_kernel_14\n");

#elif defined(STREAM_KERNEL_15)
  printf("using stream_kernel_15\n");

#elif defined(STREAM_KERNEL_16)
  printf("using stream_kernel_16\n");

#elif defined(STREAM_KERNEL_17)
  printf("using stream_kernel_17\n");

#elif defined(STREAM_KERNEL_18)
  printf("using stream_kernel_18\n");

#elif defined(STREAM_KERNEL_19)
  printf("using stream_kernel_19\n");

#elif defined(STREAM_KERNEL_20)
  printf("using stream_kernel_20\n");

#elif defined(POINTER_CHASE)
  printf("pointer chase\n");

#elif defined(C_KERNEL)
  printf("C kernel\n");

#elif defined(NAIVE_KERNEL)
  printf("using naive C stream kernel\n");

#else
  printf("no kernel\n");
  exit(1);

#endif

#ifdef STORE_DILATED_BLOCK
  printf("storing dilated A\n");
#endif

  /* Set the rand() seed. */
  if (rand_seed == 0)
  {
    srand(time(NULL));
  }

  else
  {
    srand(rand_seed);
  }

  /* Check input. */
  tree_depth = (int) ceil(log(N/(double) N_BLOCK)/log(2));
  N_padded = (int) N_BLOCK*pow(2, tree_depth);
  if (N_padded != N)
  {
    printf("N needs to be a power of 2, next larger matrix with %ix%i blocks would be %ix%i\n", N_BLOCK, N_BLOCK, N_padded, N_padded);
    N = N_padded;
  }

  if (N_only_A == 0 || N_only_A > N/N_BLOCK*N/N_BLOCK)
  {
    N_only_A = N/N_BLOCK*N/N_BLOCK;
  }

  if (N_only_B == 0 || N_only_B > N/N_BLOCK*N/N_BLOCK)
  {
    N_only_B = N/N_BLOCK*N/N_BLOCK;
  }

  if (N_only_C == 0 || N_only_C > N/N_BLOCK*N/N_BLOCK)
  {
    N_only_C = N/N_BLOCK*N/N_BLOCK;
  }

  printf("testing %ix%i matrices with a blocksize of %ix%i\n", N, N, N_BLOCK, N_BLOCK);

  /* Allocate matrices. */
#ifdef HAVE_POSIX_MEMALIGN

  if ((allocation_result = posix_memalign((void**) &A, alignment, sizeof(float)*N*N)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &B, alignment, sizeof(float)*N*N)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &C, alignment, sizeof(float)*N*N)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &D, alignment, sizeof(float)*N*N)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
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

  printf("Allocated matrices at: A = %p, B = %p, C = %p, D = %p", A, B, C, D);
#if defined(HAVE_POSIX_MEMALIGN)
  printf(", aligned on %u byte boundary", alignment);
#endif
  printf("\n");

  tree_size = sizeof(struct matrix_t);
  tree_size += sizeof(struct matrix_node_t)*pow(4, tree_depth);
  tree_size += sizeof(struct matrix_node_t*)*4*pow(4, tree_depth-1);
  tree_size += sizeof(float)*N*N;

#ifdef STORE_DILATED_BLOCK
  tree_size += sizeof(float)*N*N*4;
#endif

  if (tree_size < 1024)
  {
    printf("each matrix tree has a maximum depth of %u and a size (in case of a dense matrix) of: %llu bytes\n", tree_depth, tree_size);
  }

  else if (tree_size < 1024*1024)
  {
    printf("each matrix tree has a maximum depth of %u and a size (in case of a dense matrix) of: %1.2f kB\n", tree_depth, tree_size/1024.);
  }

  else if (tree_size < 1024*1024*1024)
  {
    printf("each matrix tree has a maximum depth of %u and a size (in case of a dense matrix) of: %1.2f MB\n", tree_depth, tree_size/1024./1024.);
  }

  else
  {
    printf("each matrix tree has a maximum depth of %u and a size (in case of a dense matrix) of: %1.2f GB\n", tree_depth, tree_size/1024./1024./1024.);
  }

  A_spamm = spamm_new(N, alignment);
  B_spamm = spamm_new(N, alignment);
  C_spamm = spamm_new(N, alignment);

  printf("Allocated spamm matrices at: A = %p, A->contiguous = %p", A_spamm, A_spamm->contiguous);
#ifdef STORE_DILATED_BLOCK
  printf(", A->contiguous_dilated = %p\n", A_spamm->contiguous_dilated);
#else
  printf("\n");
#endif

  printf("Allocated spamm matrices at: B = %p, B->contiguous = %p", B_spamm, B_spamm->contiguous);
#ifdef STORE_DILATED_BLOCK
  printf(", B->contiguous_dilated = %p\n", B_spamm->contiguous_dilated);
#else
  printf("\n");
#endif

  printf("Allocated spamm matrices at: C = %p, C->contiguous = %p", C_spamm, C_spamm->contiguous);
#ifdef STORE_DILATED_BLOCK
  printf(", C->contiguous_dilated = %p\n", C_spamm->contiguous_dilated);
#else
  printf("\n");
#endif

  /* Fill matrices with random data. In order to avoid a deterministic
   * allocation pattern, the matrix indices are reandomized as well.
   */
  index_pairs = (unsigned int (*)[2]) malloc(sizeof(unsigned int[2])*N*N);
  k = 0;
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      index_pairs[k][0] = i;
      index_pairs[k][1] = j;
      k++;
    }
  }

  if (random_elements)
  {
    printf("randomizing order of setting matrix elements\n");
    for (i = 0; i < 10*N*N; i++)
    {
      j = (int) floor(rand()/(double) RAND_MAX*N);

      k = index_pairs[i%(N*N)][0];
      index_pairs[i%(N*N)][0] = index_pairs[j][0];
      index_pairs[j][0] = k;

      k = index_pairs[i%(N*N)][1];
      index_pairs[i%(N*N)][1] = index_pairs[j][1];
      index_pairs[j][1] = k;
    }
  }

  printf("setting matrix elements\n");
  for (k = 0; k < N*N; k++)
  {
    if (random_elements)
    {
      A[get_index(N, index_pairs[k][0], index_pairs[k][1])] = rand()/(float) RAND_MAX;
      B[get_index(N, index_pairs[k][0], index_pairs[k][1])] = rand()/(float) RAND_MAX;
      C[get_index(N, index_pairs[k][0], index_pairs[k][1])] = rand()/(float) RAND_MAX;
    }

    else
    {
      A[get_index(N, index_pairs[k][0], index_pairs[k][1])] = get_index(N, index_pairs[k][0], index_pairs[k][1])+1;
      //A[get_index(N, index_pairs[k][0], index_pairs[k][1])] = (get_Morton_index(N_BLOCK, index_pairs[k][0]%N_BLOCK, index_pairs[k][1]%N_BLOCK)
      //  +interleave_2_index(index_pairs[k][0]/N_BLOCK, index_pairs[k][1]/N_BLOCK)*N_BLOCK*N_BLOCK+1)/pow(10, ceil(log(N)/log(10))+1);
      B[get_index(N, index_pairs[k][0], index_pairs[k][1])] = A[get_index(N, index_pairs[k][0], index_pairs[k][1])];
      C[get_index(N, index_pairs[k][0], index_pairs[k][1])] = A[get_index(N, index_pairs[k][0], index_pairs[k][1])];
    }

    /* Copy C matrix for verification. */
    D[get_index(N, index_pairs[k][0], index_pairs[k][1])] = C[get_index(N, index_pairs[k][0], index_pairs[k][1])];

    /* Build SpAMM tree. */
    spamm_set(A_spamm, index_pairs[k][0], index_pairs[k][1], A[get_index(N, index_pairs[k][0], index_pairs[k][1])]);
    spamm_set(B_spamm, index_pairs[k][0], index_pairs[k][1], B[get_index(N, index_pairs[k][0], index_pairs[k][1])]);
    spamm_set(C_spamm, index_pairs[k][0], index_pairs[k][1], C[get_index(N, index_pairs[k][0], index_pairs[k][1])]);
  }

  free(index_pairs);

  if (print)
  {
    printf("A =\n");
    print_matrix(N, A);
    printf("B =\n");
    print_matrix(N, B);
    printf("C =\n");
    print_matrix(N, C);

    printf("A (SpAMM) =\n");
    spamm_print(A_spamm);
    printf("B (SpAMM) =\n");
    spamm_print(B_spamm);
    printf("C (SpAMM) =\n");
    spamm_print(C_spamm);
  }

  /* Allocate stream. */
  number_stream_elements = N/N_BLOCK*N/N_BLOCK*N/N_BLOCK;
  printf("allocating stream with %llu elements\n", number_stream_elements);
  printf("1 stream element has %lu bytes\n", sizeof(struct multiply_stream_t));

  printf("mapping %u matrix blocks in stream A touching the first %lu bytes out of %lu bytes of matrix A\n",
      N_only_A, N_only_A*N_BLOCK*N_BLOCK*sizeof(float), N*N*sizeof(float));
  printf("mapping %u matrix blocks in stream B touching the first %lu bytes out of %lu bytes of matrix B\n",
      N_only_B, N_only_B*N_BLOCK*N_BLOCK*sizeof(float), N*N*sizeof(float));
  printf("mapping %u matrix blocks in stream C touching the first %lu bytes out of %lu bytes of matrix C\n",
      N_only_C, N_only_C*N_BLOCK*N_BLOCK*sizeof(float), N*N*sizeof(float));

  if (sizeof(float)*N*N < 1024)
  {
    printf("each matrix has %u blocks allocated in %lu bytes\n", N/N_BLOCK*N/N_BLOCK, sizeof(float)*N*N);
  }

  else if (sizeof(float)*N*N < 1024*1024)
  {
    printf("each matrix has %u blocks allocated in %1.2f kB\n", N/N_BLOCK*N/N_BLOCK, sizeof(float)*N*N/1024.);
  }

  else
  {
    printf("each matrix has %u blocks allocated in %1.2f MB\n", N/N_BLOCK*N/N_BLOCK, sizeof(float)*N*N/1024./1024.);
  }

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

  if (sizeof(struct multiply_stream_index_t)*number_stream_elements < 1024)
  {
    printf("multiply stream index has %llu bytes\n", sizeof(struct multiply_stream_index_t)*number_stream_elements);
  }

  else if (sizeof(struct multiply_stream_index_t)*number_stream_elements < 1024*1024)
  {
    printf("multiply stream index has %1.2f kB\n", sizeof(struct multiply_stream_index_t)*number_stream_elements/1024.);
  }

  else if (sizeof(struct multiply_stream_index_t)*number_stream_elements < 1024*1024*1024)
  {
    printf("multiply stream index has %1.2f MB\n", sizeof(struct multiply_stream_index_t)*number_stream_elements/1024./1024.);
  }

  else
  {
    printf("multiply stream index has %1.2f GB\n", sizeof(struct multiply_stream_index_t)*number_stream_elements/1024./1024./1024.);
  }

  if (sizeof(float*)*number_stream_elements < 1024)
  {
    printf("multiply stream {A,B,C} each has %llu bytes\n", sizeof(float*)*number_stream_elements);
  }

  else if (sizeof(float*)*number_stream_elements < 1024*1024)
  {
    printf("multiply stream {A,B,C} each has %1.2f kB\n", sizeof(float*)*number_stream_elements/1024.);
  }

  else if (sizeof(float*)*number_stream_elements < 1024*1024*1024)
  {
    printf("multiply stream {A,B,C} each has %1.2f MB\n", sizeof(float*)*number_stream_elements/1024./1024.);
  }

  else
  {
    printf("multiply stream {A,B,C} each has %1.2f GB\n", sizeof(float*)*number_stream_elements/1024./1024./1024.);
  }

#ifdef HAVE_POSIX_MEMALIGN
  if ((allocation_result = posix_memalign((void**) &multiply_stream, alignment, sizeof(struct multiply_stream_t)*number_stream_elements)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &A_stream, alignment, sizeof(float*)*number_stream_elements)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &B_stream, alignment, sizeof(float*)*number_stream_elements)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &C_stream, alignment, sizeof(float*)*number_stream_elements)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if ((allocation_result = posix_memalign((void**) &unique_blocks, alignment, sizeof(void*)*(N_only_A+N_only_B+N_only_C))) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }

  if (histogram)
  {
    if ((allocation_result = posix_memalign((void**) &histogram_distances, alignment, sizeof(unsigned long long)*(N_only_A+N_only_B+N_only_C)*(N_only_A+N_only_B+N_only_C))) != 0)
    {
      switch (allocation_result)
      {
        case EINVAL:
          printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
          break;

        case ENOMEM:
          printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
          break;

        default:
          printf("Unknown error code\n");
          break;
      }
      exit(1);
    }

    if ((allocation_result = posix_memalign((void**) &histogram_counts, alignment, sizeof(unsigned int)*(N_only_A+N_only_B+N_only_C)*(N_only_A+N_only_B+N_only_C))) != 0)
    {
      switch (allocation_result)
      {
        case EINVAL:
          printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
          break;

        case ENOMEM:
          printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
          break;

        default:
          printf("Unknown error code\n");
          break;
      }
      exit(1);
    }
  }

  if ((allocation_result = posix_memalign((void**) &multiply_stream_index, alignment, sizeof(struct multiply_stream_t)*number_stream_elements)) != 0)
  {
    switch (allocation_result)
    {
      case EINVAL:
        printf("The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n");
        break;

      case ENOMEM:
        printf("[line %u] There was insufficient memory to fulfill the allocation request.\n", __LINE__);
        break;

      default:
        printf("Unknown error code\n");
        break;
    }
    exit(1);
  }
#else
  if ((multiply_stream = (struct multiply_stream_t*) malloc(sizeof(struct multiply_stream_t)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream\n");
    exit(1);
  }

  if ((A_stream = (float**) malloc(sizeof(float*)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream A\n");
    exit(1);
  }

  if ((B_stream = (float**) malloc(sizeof(float*)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream B\n");
    exit(1);
  }

  if ((C_stream = (float**) malloc(sizeof(float*)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream C\n");
    exit(1);
  }

  if ((unique_blocks = (void**) malloc(sizeof(void*)*(N_only_A+N_only_B+N_only_C))) == NULL)
  {
    printf("error allocating list of unique blocks\n");
    exit(1);
  }

  if (histogram)
  {
    if ((histogram_distances = (unsigned long long*) malloc(sizeof(unsigned long long)*(N_only_A+N_only_B+N_only_C)*(N_only_A+N_only_B+N_only_C))) == NULL)
    {
      printf("error allocating list of histogram distances\n");
      exit(1);
    }

    if ((histogram_counts = (unsigned int*) malloc(sizeof(unsigned int)*(N_only_A+N_only_B+N_only_C)*(N_only_A+N_only_B+N_only_C))) == NULL)
    {
      printf("error allocating list of histogram counts\n");
      exit(1);
    }
  }

  if ((multiply_stream_index = (struct multiply_stream_index_t*) malloc(sizeof(struct multiply_stream_index_t)*number_stream_elements)) == NULL)
  {
    printf("error allocating multiply stream index\n");
    exit(1);
  }
#endif

  printf("allocated multiply stream at %p\n", multiply_stream);
  printf("allocated multiply stream A at %p\n", A_stream);
  printf("allocated multiply stream B at %p\n", B_stream);
  printf("allocated multiply stream C at %p\n", C_stream);
  printf("allocated multiply stream index at %p\n", multiply_stream_index);

  stream_index = 0;
  unique_block_index = 0;
  /* In order to avoid cache aliasing effects between blocks, we increment the
   * 3 stream indices starting at different values.
   */
  stream_index_A = 0;
  stream_index_B = 1;
  stream_index_C = 2;

#ifdef DENSE_MULTIPLY
  for (i = 0; i < N; i += N_BLOCK) {
    for (j = 0; j < N; j += N_BLOCK) {
      for (k = 0; k < N; k += N_BLOCK)
      {
        if (stream_index < N_only_A)
        {
          multiply_stream[stream_index].A_block = &A[get_blocked_index(N, i, k)];
          A_stream[stream_index] = &A[get_blocked_index(N, i, k)];
          multiply_stream_index[stream_index].A_index = get_blocked_index(N, i, k)/N_BLOCK/N_BLOCK;
          unique_blocks[unique_block_index++] = &A[get_blocked_index(N, i, k)];
        }

        else
        {
          multiply_stream[stream_index].A_block = multiply_stream[stream_index_A%N_only_A].A_block;
          A_stream[stream_index] = A_stream[stream_index_A%N_only_A];
          multiply_stream_index[stream_index].A_index = multiply_stream_index[stream_index_A%N_only_A].A_index;
        }

        if (stream_index < N_only_B)
        {
          multiply_stream[stream_index].B_block = &B[get_blocked_index(N, k, j)];
          B_stream[stream_index] = &B[get_blocked_index(N, k, j)];
          multiply_stream_index[stream_index].B_index = get_blocked_index(N, k, j)/N_BLOCK/N_BLOCK;
          unique_blocks[unique_block_index++] = &B[get_blocked_index(N, k, j)];
        }

        else
        {
          multiply_stream[stream_index].B_block = multiply_stream[stream_index_B%N_only_B].B_block;
          B_stream[stream_index] = B_stream[stream_index_B%N_only_B];
          multiply_stream_index[stream_index].B_index = multiply_stream_index[stream_index_B%N_only_B].B_index;
        }

        if (stream_index < N_only_C)
        {
          multiply_stream[stream_index].C_block = &C[get_blocked_index(N, i, j)];
          C_stream[stream_index] = &C[get_blocked_index(N, i, j)];
          multiply_stream_index[stream_index].C_index = get_blocked_index(N, i, j)/N_BLOCK/N_BLOCK;
          unique_blocks[unique_block_index++] = &C[get_blocked_index(N, i, j)];
        }

        else
        {
          multiply_stream[stream_index].C_block = multiply_stream[stream_index_C%N_only_C].C_block;
          C_stream[stream_index] = C_stream[stream_index_C%N_only_C];
          multiply_stream_index[stream_index].C_index = multiply_stream_index[stream_index_C%N_only_C].C_index;
        }

        multiply_stream_index[stream_index].index = interleave_3_index(multiply_stream_index[stream_index].A_index,
            multiply_stream_index[stream_index].B_index,
            multiply_stream_index[stream_index].C_index);

        stream_index++;

        /* Avoid cache aliasing effects. */
        if (stream_index%N_only_A == 0) { stream_index_A++; }
        if (stream_index%N_only_B == 0) { stream_index_B++; }
        if (stream_index%N_only_C == 0) { stream_index_C++; }

        stream_index_A++;
        stream_index_B++;
        stream_index_C++;
      }
    }
  }
#else
  spamm_multiply(A_spamm, B_spamm, C_spamm, &stream_index, multiply_stream, multiply_stream_index, A_stream, B_stream, C_stream);
#endif

  printf("created %llu multiply stream elements\n", stream_index);

  if (print)
  {
    printf("multiply_stream_index = ");
    for (i = 0; i < number_stream_elements; i++)
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
      for (i = 0; i < number_stream_elements; i++)
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

  /* Calculate histogram of address distances of all blocks on stream. */
  if (histogram)
  {
    /* Sort and uniqify list of unique blocks (to really make it a list of
     * unique blocks).
     */
    for (i = 0; i < N_only_A+N_only_B+N_only_C-1; i++) {
      for (j = i; j < N_only_A+N_only_B+N_only_C; j++)
      {
        if (unique_blocks[i] > unique_blocks[j])
        {
          temp_pointer = unique_blocks[i];
          unique_blocks[i] = unique_blocks[j];
          unique_blocks[j] = temp_pointer;
        }
      }
    }

    k = 0;
    for (i = 0; i < N_only_A+N_only_B+N_only_C; i++) {
      for (j = i+1; j < N_only_A+N_only_B+N_only_C; j++)
      {
        if (unique_blocks[i] != unique_blocks[j])
        {
          histogram_distances[k++] = (intptr_t) unique_blocks[j] - (intptr_t) unique_blocks[i];
        }
      }

      /* Skip over identical blocks. */
      while (i < N_only_A+N_only_B+N_only_C-1 && unique_blocks[i] == unique_blocks[i+1]) { i++; }
    }

    /* Count. */
    for (i = 0; i < k-1; i++) {
      for (j = i+1; j < k; j++)
      {
        if (histogram_distances[i] > histogram_distances[j])
        {
          temp_unsigned_int = histogram_distances[i];
          histogram_distances[i] = histogram_distances[j];
          histogram_distances[j] = temp_unsigned_int;
        }
      }
    }

    //printf("unique block address distance histogram (only multiples of 4096)\n");
    //printf("distance in bytes, count\n");
    temp_unsigned_int = 0;
    for (i = 0; i < k; i++) {
      histogram_counts[i] = 1;
      j = i;
      while (i+1 < k && histogram_distances[i] == histogram_distances[i+1])
      {
        histogram_counts[j]++;
        i++;
      }
      if (histogram_distances[j]%4096 == 0)
      {
        temp_unsigned_int += histogram_counts[j];
        //printf("%llu %u\n", histogram_distances[j], histogram_counts[j]);
      }
    }
    printf("unique block address distance multiples of 4096 = %u\n", temp_unsigned_int);
  }

  /* Apply beta to C. */
#if ! defined(EXTERNAL_BLAS)
  for (i = 0; i < N; i++) {
    for (j = 0; j < N; j++)
    {
      C[i*N+j] *= beta;
      spamm_set(C_spamm, i, j, beta*spamm_get(C_spamm, i, j));
      if (verify) { D[i*N+j] *= beta; }
    }
  }
  printf("applied beta\n");

  if (print)
  {
    printf("C =\n");
    print_matrix(N, C);

    printf("C (SpAMM) =\n");
    spamm_print(C_spamm);
  }
#endif

  printf("looping over multiply (loops = %llu)\n", loops);

  gettimeofday(&start, NULL);
  getrusage(RUSAGE_SELF, &rusage_start);

#ifdef HAVE_PAPI
  papi_errorcode = PAPI_num_events(papi_events);
  if (papi_errorcode == 0)
  {
    printf("[PAPI] no counters specified\n");
  }

  else if (papi_errorcode < 0)
  {
    print_papi_error("failed to count number of counters in eventset", papi_errorcode);
  }

  else
  {
    if ((papi_errorcode = PAPI_start(papi_events))) { print_papi_error("failed start", papi_errorcode); }
  }
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
    stream_multiply(number_stream_elements, alpha, multiply_stream, A_stream, B_stream, C_stream);

#endif
  }
#endif

#ifdef HAVE_PAPI
  if (PAPI_num_events(papi_events) > 0)
  {
    if (PAPI_stop(papi_events, papi_values) != PAPI_OK) { printf("failed to stop\n"); exit(1); }
    papi_value_index = 0;

    if (load_TOT_INS)
    {
      printf("[PAPI] total instructions = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_TOT_CYC)
    {
      printf("[PAPI] total cycles = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_RES_STL)
    {
      printf("[PAPI] cycles stalled = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L1_ICM)
    {
      printf("[PAPI] instruction L1 misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L1_DCM)
    {
      printf("[PAPI] data L1 misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L1_LDM)
    {
      printf("[PAPI] data L1 load misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L1_STM)
    {
      printf("[PAPI] data L1 store misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L1_DCH)
    {
      printf("[PAPI] data L1 hit = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L1_DCA)
    {
      printf("[PAPI] data L1 access = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L2_ICM)
    {
      printf("[PAPI] instruction L2 misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L2_DCM)
    {
      printf("[PAPI] data L2 misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L2_LDM)
    {
      printf("[PAPI] data L2 load misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L2_STM)
    {
      printf("[PAPI] data L2 store misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L2_DCH)
    {
      printf("[PAPI] data L2 hit = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_L2_DCA)
    {
      printf("[PAPI] data L2 access = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_TLB_IM)
    {
      printf("[PAPI] total instruction TLB misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_TLB_DM)
    {
      printf("[PAPI] total data TLB misses = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }

    if (load_TLB_SD)
    {
      printf("[PAPI] total data TLB shootdowns = %lli per stream element = %1.2f\n",
          papi_values[papi_value_index], papi_values[papi_value_index]/(double) loops/(double) number_stream_elements);
      papi_value_index++;
    }
  }

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
    printf("result: C =\n");
    print_matrix(N, C);
    printf("result: C_spamm =\n");
    spamm_print(C_spamm);
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

#ifdef DENSE_MULTIPLY
    max_diff = compare_matrix(N, C, D);
#else
    max_diff = spamm_compare_matrix(C_spamm, D);
#endif
    printf("max diff = %e\n", max_diff);
  }

  /* Free memory. */
  free(A);
  free(B);
  free(C);
  free(D);
  spamm_free(&A_spamm);
  spamm_free(&B_spamm);
  spamm_free(&C_spamm);
  free(multiply_stream);
  free(multiply_stream_index);
}
