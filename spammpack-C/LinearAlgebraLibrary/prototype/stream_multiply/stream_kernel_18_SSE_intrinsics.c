#include <stdlib.h>
#include <stdio.h>
#include <xmmintrin.h>

#define A_OFFSET_11 (0*4+0)*64
#define A_OFFSET_12 (0*4+1)*64
#define A_OFFSET_13 (0*4+2)*64
#define A_OFFSET_14 (0*4+3)*64
#define A_OFFSET_21 (1*4+0)*64
#define A_OFFSET_22 (1*4+1)*64
#define A_OFFSET_23 (1*4+2)*64
#define A_OFFSET_24 (1*4+3)*64
#define A_OFFSET_31 (2*4+0)*64
#define A_OFFSET_32 (2*4+1)*64
#define A_OFFSET_33 (2*4+2)*64
#define A_OFFSET_34 (2*4+3)*64
#define A_OFFSET_41 (3*4+0)*64
#define A_OFFSET_42 (3*4+1)*64
#define A_OFFSET_43 (3*4+2)*64
#define A_OFFSET_44 (3*4+3)*64

#define B_OFFSET_11 (0*0+0)*16
#define B_OFFSET_12 (0*0+1)*16
#define B_OFFSET_13 (0*0+2)*16
#define B_OFFSET_14 (0*0+3)*16
#define B_OFFSET_21 (1*0+0)*16
#define B_OFFSET_22 (1*0+1)*16
#define B_OFFSET_23 (1*0+2)*16
#define B_OFFSET_24 (1*0+3)*16
#define B_OFFSET_31 (2*0+0)*16
#define B_OFFSET_32 (2*0+1)*16
#define B_OFFSET_33 (2*0+2)*16
#define B_OFFSET_34 (2*0+3)*16
#define B_OFFSET_41 (3*0+0)*16
#define B_OFFSET_42 (3*0+1)*16
#define B_OFFSET_43 (3*0+2)*16
#define B_OFFSET_44 (3*0+3)*16

#define C_OFFSET_11 (0*0+0)*16
#define C_OFFSET_12 (0*0+1)*16
#define C_OFFSET_13 (0*0+2)*16
#define C_OFFSET_14 (0*0+3)*16
#define C_OFFSET_21 (1*0+0)*16
#define C_OFFSET_22 (1*0+1)*16
#define C_OFFSET_23 (1*0+2)*16
#define C_OFFSET_24 (1*0+3)*16
#define C_OFFSET_31 (2*0+0)*16
#define C_OFFSET_32 (2*0+1)*16
#define C_OFFSET_33 (2*0+2)*16
#define C_OFFSET_34 (2*0+3)*16
#define C_OFFSET_41 (3*0+0)*16
#define C_OFFSET_42 (3*0+1)*16
#define C_OFFSET_43 (3*0+2)*16
#define C_OFFSET_44 (3*0+3)*16

struct multiply_stream_t
{
  float *A_block;
  float *B_block;
  float *C_block;
  float  norm[32];
};

void
stream_kernel_18_SSE_intrinsics (const unsigned int number_stream_elements,
    float alpha,
    float tolerance,
    struct multiply_stream_t *multiply_stream)
{
  short int i;
  unsigned int stream_index;
  unsigned int max_stream_index;

  float *restrict A;
  float *restrict B;
  float *restrict C;

  float *restrict norm;

  __m128 alpha_row;

  __m128 A_element;
  __m128 B_row;
  __m128 C_row[4];


  /* Divide number of stream elements by 64 to simulate stride of 64. */
  max_stream_index = number_stream_elements/64;

  alpha_row = _mm_set1_ps(alpha);

  for (stream_index = 0; stream_index < max_stream_index; stream_index++)
  {
    /* Load pointers to matrix data blocks. */
    A = multiply_stream[stream_index].A_block;
    B = multiply_stream[stream_index].B_block;
    C = multiply_stream[stream_index].C_block;
    norm = multiply_stream[stream_index].norm;

    if (norm[0]*norm[4] >= tolerance && norm[1]*norm[6] >= tolerance)
    {
      /* A(1,1)*B(1,1) = C(1,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      /* A(1,2)*B(2,1) = C(1,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      for (i = 0; i < 4; i++)
      {
        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
      }
    }

    else if (norm[0]*norm[4] >= tolerance)
    {
      /* A(1,1)*B(1,1) = C(1,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
      }
    }

    else if (norm[1]*norm[6] >= tolerance)
    {
      /* A(1,2)*B(2,1) = C(1,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_11]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_11], C_row[i]);
      }
    }

    if (norm[2]*norm[4] >= tolerance && norm[3]*norm[6] >= tolerance)
    {
      /* A(2,1)*B(1,1) = C(2,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      /* A(2,2)*B(2,1) = C(2,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      for (i = 0; i < 4; i++)
      {
        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
      }
    }

    else if (norm[2]*norm[4] >= tolerance)
    {
      /* A(2,1)*B(1,1) = C(2,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_11]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_11]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
      }
    }

    else if (norm[3]*norm[6] >= tolerance)
    {
      /* A(2,2)*B(2,1) = C(2,1). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_21]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_21]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_21]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_21], C_row[i]);
      }
    }

    if (norm[0]*norm[5] >= tolerance && norm[1]*norm[7] >= tolerance)
    {
      /* A(1,1)*B(1,2) = C(1,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      /* A(1,2)*B(2,2) = C(1,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      for (i = 0; i < 4; i++)
      {
        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
      }
    }

    else if (norm[0]*norm[5] >= tolerance)
    {
      /* A(1,1)*B(1,2) = C(1,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_11]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
      }
    }

    else if (norm[1]*norm[7] >= tolerance)
    {
      /* A(1,2)*B(2,2) = C(1,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_12]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_12]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_12], C_row[i]);
      }
    }

    if (norm[2]*norm[5] >= tolerance && norm[3]*norm[7] >= tolerance)
    {
      /* A(2,1)*B(1,2) = C(2,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      /* A(2,2)*B(2,2) = C(2,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);
      }

      for (i = 0; i < 4; i++)
      {
        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
      }
    }

    else if (norm[2]*norm[5] >= tolerance)
    {
      /* A(2,1)*B(1,2) = C(2,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_12]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_21]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_12]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
      }
    }

    else if (norm[3]*norm[7] >= tolerance)
    {
      /* A(2,2)*B(2,2) = C(2,2). */
      for (i = 0; i < 4; i++)
      {
        A_element = _mm_load_ps(&A[(i*4+0)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[0*4+B_OFFSET_22]);
        C_row[i] = _mm_mul_ps(A_element, B_row);

        A_element = _mm_load_ps(&A[(i*4+1)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[1*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+2)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[2*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        A_element = _mm_load_ps(&A[(i*4+3)*4+A_OFFSET_22]);
        B_row = _mm_load_ps(&B[3*4+B_OFFSET_22]);
        C_row[i] = _mm_add_ps(_mm_mul_ps(A_element, B_row), C_row[i]);

        C_row[i] = _mm_mul_ps(alpha_row, C_row[i]);
        C_row[i] = _mm_add_ps(_mm_load_ps(&C[i*4+C_OFFSET_22]), C_row[i]);
        _mm_store_ps(&C[i*4+C_OFFSET_22], C_row[i]);
      }
    }
  }
}
