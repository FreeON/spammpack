#include <stdlib.h>
#include <xmmintrin.h>

#define A_OFFSET_11 0
#define A_OFFSET_12 64
#define A_OFFSET_21 128
#define A_OFFSET_22 192

#define B_OFFSET_11 0
#define B_OFFSET_12 16
#define B_OFFSET_21 32
#define B_OFFSET_22 64

#define C_OFFSET_11 0
#define C_OFFSET_12 16
#define C_OFFSET_21 32
#define C_OFFSET_22 64

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
  char  mask[8];
};

void
stream_kernel_15_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    struct multiply_stream_t *multiply_stream)
{
  short int i;
  unsigned int stream_index;
  unsigned int max_stream_index;

  float *restrict A;
  float *restrict B;
  float *restrict C;

  char *restrict mask;

  __m128 alpha_row;

  __m128 C_row_1;
  __m128 C_row_2;
  __m128 C_row_3;
  __m128 C_row_4;

  /* Divide number of stream elements by 8 to simulate stride of 8. */
  max_stream_index = number_stream_elements/8;

  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
  {
    /* Load pointers to matrix data blocks. */
    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;
    mask = multiply_stream[stream_index].mask;

    C_row_1 = _mm_setzero_ps();
    C_row_2 = _mm_setzero_ps();
    C_row_3 = _mm_setzero_ps();
    C_row_4 = _mm_setzero_ps();

    if (mask[0]) /* A(1,1)*B(1,1) = C(1,1). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_11]), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_11]), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_11), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_11), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_11]), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_11]), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_11]), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_11]), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_4))));
    }

    if (mask[1]) /* A(1,2)*B(2,1) = C(1,1). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_12]), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_12]), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_12), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_12), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_12]), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_12]), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_12]), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_12]), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_4))));
    }

    if (mask[0] | mask[1])
    {
      C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
      C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
      C_row_3 = _mm_mul_ps(alpha_row, C_row_3);
      C_row_4 = _mm_mul_ps(alpha_row, C_row_4);

      C_row_1 = _mm_add_ps(_mm_load_ps(&C[0*4+C_OFFSET_11]), C_row_1);
      C_row_2 = _mm_add_ps(_mm_load_ps(&C[1*4+C_OFFSET_11]), C_row_2);
      C_row_3 = _mm_add_ps(_mm_load_ps(&C[2*4+C_OFFSET_11]), C_row_3);
      C_row_4 = _mm_add_ps(_mm_load_ps(&C[3*4+C_OFFSET_11]), C_row_4);

      _mm_store_ps(&C[0*4+C_OFFSET_11], C_row_1);
      _mm_store_ps(&C[1*4+C_OFFSET_11], C_row_2);
      _mm_store_ps(&C[2*4+C_OFFSET_11], C_row_3);
      _mm_store_ps(&C[3*4+C_OFFSET_11], C_row_4);
    }

    C_row_1 = _mm_setzero_ps();
    C_row_2 = _mm_setzero_ps();
    C_row_3 = _mm_setzero_ps();
    C_row_4 = _mm_setzero_ps();

    if (mask[2]) /* A(2,1)*B(1,1) = C(2,1). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_21]), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_21]), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_21), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_21), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_21]), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_21]), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_21]), _mm_load_ps(&B[3*4+B_OFFSET_11])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_21]), _mm_load_ps(&B[2*4+B_OFFSET_11])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_11])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_11])),
                C_row_4))));
    }

    if (mask[3]) /* A(2,2)*B(2,1) = C(2,1). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_22]), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_22]), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_22), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_22), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_22]), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_22]), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_22]), _mm_load_ps(&B[3*4+B_OFFSET_21])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_22]), _mm_load_ps(&B[2*4+B_OFFSET_21])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_21])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_21])),
                C_row_4))));
    }

    if (mask[2] | mask[3])
    {
      C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
      C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
      C_row_3 = _mm_mul_ps(alpha_row, C_row_3);
      C_row_4 = _mm_mul_ps(alpha_row, C_row_4);

      C_row_1 = _mm_add_ps(_mm_load_ps(&C[0*4+C_OFFSET_21]), C_row_1);
      C_row_2 = _mm_add_ps(_mm_load_ps(&C[1*4+C_OFFSET_21]), C_row_2);
      C_row_3 = _mm_add_ps(_mm_load_ps(&C[2*4+C_OFFSET_21]), C_row_3);
      C_row_4 = _mm_add_ps(_mm_load_ps(&C[3*4+C_OFFSET_21]), C_row_4);

      _mm_store_ps(&C[0*4+C_OFFSET_21], C_row_1);
      _mm_store_ps(&C[1*4+C_OFFSET_21], C_row_2);
      _mm_store_ps(&C[2*4+C_OFFSET_21], C_row_3);
      _mm_store_ps(&C[3*4+C_OFFSET_21], C_row_4);
    }

    C_row_1 = _mm_setzero_ps();
    C_row_2 = _mm_setzero_ps();
    C_row_3 = _mm_setzero_ps();
    C_row_4 = _mm_setzero_ps();

    if (mask[4]) /* A(1,1)*B(1,2) = C(1,2). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_11]), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_11]), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_11), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_11), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_11]), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_11]), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_11]), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_11]), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_11]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_11]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_4))));
    }

    if (mask[5]) /* A(1,2)*B(2,2) = C(1,2). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_12]), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_12]), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_12), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_12), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_12]), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_12]), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_12]), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_12]), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_12]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_12]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_4))));
    }

    if (mask[4] | mask[5])
    {
      C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
      C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
      C_row_3 = _mm_mul_ps(alpha_row, C_row_3);
      C_row_4 = _mm_mul_ps(alpha_row, C_row_4);

      C_row_1 = _mm_add_ps(_mm_load_ps(&C[0*4+C_OFFSET_12]), C_row_1);
      C_row_2 = _mm_add_ps(_mm_load_ps(&C[1*4+C_OFFSET_12]), C_row_2);
      C_row_3 = _mm_add_ps(_mm_load_ps(&C[2*4+C_OFFSET_12]), C_row_3);
      C_row_4 = _mm_add_ps(_mm_load_ps(&C[3*4+C_OFFSET_12]), C_row_4);

      _mm_store_ps(&C[0*4+C_OFFSET_12], C_row_1);
      _mm_store_ps(&C[1*4+C_OFFSET_12], C_row_2);
      _mm_store_ps(&C[2*4+C_OFFSET_12], C_row_3);
      _mm_store_ps(&C[3*4+C_OFFSET_12], C_row_4);
    }

    C_row_1 = _mm_setzero_ps();
    C_row_2 = _mm_setzero_ps();
    C_row_3 = _mm_setzero_ps();
    C_row_4 = _mm_setzero_ps();

    if (mask[6]) /* A(2,1)*B(1,2) = C(2,2). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_21]), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_21]), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_21), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_21), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_21]), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_21]), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_21]), _mm_load_ps(&B[3*4+B_OFFSET_12])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_21]), _mm_load_ps(&B[2*4+B_OFFSET_12])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_21]), _mm_load_ps(&B[1*4+B_OFFSET_12])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_21]), _mm_load_ps(&B[0*4+B_OFFSET_12])),
                C_row_4))));
    }

    if (mask[7]) /* A(2,2)*B(2,2) = C(2,2). */
    {
      C_row_1 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+3)*4+A_OFFSET_22]), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+2)*4+A_OFFSET_22]), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(0*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_1))));

      C_row_2 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+3)*4]+A_OFFSET_22), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+2)*4]+A_OFFSET_22), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(1*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_2))));

      C_row_3 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+3)*4+A_OFFSET_22]), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+2)*4+A_OFFSET_22]), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(2*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_3))));

      C_row_4 = _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+3)*4+A_OFFSET_22]), _mm_load_ps(&B[3*4+B_OFFSET_22])),
          _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+2)*4+A_OFFSET_22]), _mm_load_ps(&B[2*4+B_OFFSET_22])),
            _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+1)*4+A_OFFSET_22]), _mm_load_ps(&B[1*4+B_OFFSET_22])),
              _mm_add_ps(_mm_mul_ps(_mm_load_ps(&A[(3*4+0)*4+A_OFFSET_22]), _mm_load_ps(&B[0*4+B_OFFSET_22])),
                C_row_4))));
    }

    if (mask[6] | mask[7])
    {
      C_row_1 = _mm_mul_ps(alpha_row, C_row_1);
      C_row_2 = _mm_mul_ps(alpha_row, C_row_2);
      C_row_3 = _mm_mul_ps(alpha_row, C_row_3);
      C_row_4 = _mm_mul_ps(alpha_row, C_row_4);

      C_row_1 = _mm_add_ps(_mm_load_ps(&C[0*4+C_OFFSET_22]), C_row_1);
      C_row_2 = _mm_add_ps(_mm_load_ps(&C[1*4+C_OFFSET_22]), C_row_2);
      C_row_3 = _mm_add_ps(_mm_load_ps(&C[2*4+C_OFFSET_22]), C_row_3);
      C_row_4 = _mm_add_ps(_mm_load_ps(&C[3*4+C_OFFSET_22]), C_row_4);

      _mm_store_ps(&C[0*4+C_OFFSET_22], C_row_1);
      _mm_store_ps(&C[1*4+C_OFFSET_22], C_row_2);
      _mm_store_ps(&C[2*4+C_OFFSET_22], C_row_3);
      _mm_store_ps(&C[3*4+C_OFFSET_22], C_row_4);
    }
  }
}
